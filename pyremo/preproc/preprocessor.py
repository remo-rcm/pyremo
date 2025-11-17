import os
import tempfile

import dask
import xarray as xr

import pyremo as pr
from warnings import warn
from ..remo_ds import update_meta_info, parse_dates
from .era5 import era5_gfile_from_dkrz
from .utils import datelist, ensure_dir, write_forcing_file
from .remapping import remap, remap_remo
from .cf import get_gcm_dataset, get_gcm_gfile
from ..tables import domains


dkrz_template = {
    "path_template": "/pool/data/ERA5/{era_id}/{level_type}/{dataType}/{frequency}/{code:03d}",
    "file_template": "{era_id}{level_type}{typeid}_{frequency}_{date}_{code:03d}.grb",
}

cloud_urls = {
    "gc": "https://storage.googleapis.com/cmip6/pangeo-cmip6.json",
    "aws": "https://cmip6-pds.s3.amazonaws.com/pangeo-cmip6.json",
}


def open_zstore(zstore):
    """
    Open a Zarr store.

    Parameters
    ----------
    zstore : str
        The path to the Zarr store.

    Returns
    -------
    xarray.Dataset
        The opened dataset.
    """
    import fsspec

    # logger.debug(f"opening: {zstore}")
    return xr.open_zarr(store=fsspec.get_mapper(zstore, anon=True), consolidated=True)


def get_dsets(df):
    """
    Get datasets from a DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing the dataset information.

    Returns
    -------
    dict
        Dictionary of datasets.
    """
    dsets = {row.variable_id: open_zstore(row.zstore) for row in df.itertuples()}
    return dsets


def merge_dsets(dsets, ref="ta"):
    """
    Merge datasets with a reference dataset.

    Parameters
    ----------
    dsets : dict
        Dictionary of datasets to merge.
    ref : str, optional
        Reference variable ID, by default "ta".

    Returns
    -------
    list
        List of merged datasets.
    """
    ref_ds = dsets[ref]
    dsets_list = [ref_ds.convert_calendar(ref_ds.time.dt.calendar, use_cftime=True)]
    for variable_id, ds in dsets.items():
        if "time" in ds.coords:
            ds = ds.convert_calendar(ds.time.dt.calendar, use_cftime=True)
        if variable_id != ref:
            dsets_list.append(ds[variable_id])
    return dsets_list


def get_pangeo_catalog(
    url,
    source_id,
    experiment_id="historical",
    member_id="r1i1p1f1",
    activity_id=None,
    **kwargs,
):
    """
    Get Zarr datasets from the cloud.

    Parameters
    ----------
    url : str
        URL of the cloud storage.
    source_id : str
        Source ID of the dataset.
    experiment_id : str, optional
        Experiment ID, by default "historical".
    member_id : str, optional
        Member ID, by default "r1i1p1f1".
    activity_id : str, optional
        Activity ID, by default None.

    Returns
    -------
    intake.catalog
        Intake catalog of the searched datasets.
    """
    import intake

    catalog = intake.open_esm_datastore(url)
    if activity_id is None:
        activity_id = "CMIP" if experiment_id == "historical" else "ScenarioMIP"

    return catalog.search(
        activity_id=activity_id,
        experiment_id=experiment_id,
        table_id=["6hrLev", "fx", "Oday"],
        variable_id=["ta", "ua", "va", "hus", "sftlf", "orog", "tos"],
        member_id=member_id,
        source_id=source_id,
        **kwargs,
    )


def get_filename(date, expid="000000", template=None):
    """
    Get the filename for a given date and experiment ID.

    Parameters
    ----------
    date : datetime-like
        The date for the filename.
    expid : str, optional
        Experiment ID, by default "000000".
    template : str, optional
        Filename template, by default None.

    Returns
    -------
    str
        Formatted filename.
    """
    if template is None:
        template = "a{expid}a{date:%Y%m%d%H}.nc"
    return template.format(expid=expid, date=date)


def get_outpath(path, time):
    """
    Get the output path for a given time.

    Parameters
    ----------
    path : str
        Path template.
    time : datetime-like
        Time for the output path.

    Returns
    -------
    str
        Formatted output path.
    """
    outpath = path.format(date=time)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    return outpath


def prepare_surflib(surflib):
    """
    Prepare the surface library.

    Parameters
    ----------
    surflib : str
        Path to the surface library.

    Returns
    -------
    xarray.Dataset
        Prepared surface library dataset.
    """
    return update_meta_info(xr.open_dataset(surflib).squeeze(drop=True).load())


class Preprocessor:
    """Generic preprocessing pipeline for preparing forcing input data.

    This base class encapsulates common logic for:

    * Managing a temporary scratch directory.
    * Loading / preparing a surface library (``surflib``).
    * Initializing domain metadata (either from a registry or inferred).
    * Remapping / transforming raw input datasets into model forcing.
    * Writing forcing files for a sequence of timesteps.

    Subclasses override :meth:`get_input_dataset` and optionally ``self.remap``
    to adapt to different upstream data sources (CF-compliant, REMO, ERA5, cloud).

    Parameters
    ----------
    expid : str
        Experiment identifier used in output naming.
    surflib : str
        Path to surface library NetCDF file.
    domain : dict or str, optional
        Domain metadata dictionary or a registered domain id. If ``None`` the
        domain is inferred from ``surflib`` with a 1-grid-cell interior crop.
    vc : str or dict, optional
        Vertical coordinate table key (looked up in ``pr.vc.tables``) or an
        explicit vertical coordinate mapping. Defaults to ``"vc_49lev"``.
    outpath : str, optional
        Output path template (e.g. ``"/data/run/{date:%Y%m%d}"``). If set, it
        is used by :meth:`run` when writing forcing files.
    scratch : str, optional
        Parent directory for an auto-created temporary working directory. If
        ``None`` uses the ``SCRATCH`` environment variable.

    Attributes
    ----------
    expid : str
        Experiment identifier.
    surflib : xarray.Dataset
        Prepared surface library dataset.
    domain_info : dict
        Domain metadata used during remapping.
    vc : dict
        Vertical coordinate mapping / table.
    outpath : str or None
        Output path template.
    scratch : tempfile.TemporaryDirectory
        Temporary working directory context.
    remap : callable
        Function applied to raw datasets producing forcing variables.
    """

    def __init__(
        self, expid, surflib, domain=None, vc=None, outpath=None, scratch=None
    ):
        if scratch is None:
            scratch = os.environ["SCRATCH"]
        else:
            os.makedirs(scratch, exist_ok=True)
        self.scratch = tempfile.TemporaryDirectory(dir=scratch)
        print(f"Preprocessor scratch: {self.scratch.name}")
        # logger.debug(f"scratch: {self.scratch.name}")
        if vc is None:
            vc = "vc_49lev"
        if isinstance(vc, str):
            self.vc = pr.vc.tables[vc]
        else:
            self.vc = vc
        self.expid = expid
        self.surflib = prepare_surflib(surflib)
        self.outpath = outpath
        self.remap = remap
        self._init_domain_info(domain)

    def _clean(self):
        """Remove the temporary scratch directory.

        Called automatically at the end of :meth:`run` when results have been
        computed and written to disk.
        """
        # logger.debug(f"cleaning up: {self.scratch.name}")
        self.scratch.cleanup()

    def _init_domain_info(self, domain=None):
        """Initialize domain metadata used for horizontal remapping.

        Domain information can be provided explicitly as a dict, referenced by
        a registered domain id, or inferred from the ``surflib`` dataset. When
        inferring, one grid cell is cropped from each outer boundary and the
        lower-left corner adjusted accordingly.

        Parameters
        ----------
        domain : dict or str or None, optional
            Domain metadata or domain id. If ``None`` infer from ``surflib``.
        """
        if isinstance(domain, dict):
            self.domain_info = domain
        elif isinstance(domain, str) and domain in domains.table.index:
            self.domain_info = pr.domain_info(domain)
        else:
            warn("domain_id not registered, taking grid from surflib...")
            domain_info = self.surflib.cx.info()
            domain_info["nlon"] -= 2
            domain_info["nlat"] -= 2
            domain_info["ll_lon"] += domain_info["dlon"]
            domain_info["ll_lat"] += domain_info["dlat"]
            domain_info["domain_id"] = domain or "custom"
            self.domain_info = domain_info

    def write(self, ds, outpath):
        """Write a forcing dataset to NetCDF.

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset containing all required forcing variables.
        outpath : str
            Directory path where the file will be created.

        Returns
        -------
        str
            Absolute path to the written forcing file.
        """
        return write_forcing_file(ds, path=outpath, expid=self.expid)

    @dask.delayed
    def preprocess(self, date=None, ds=None, outpath=None, initial=None, **kwargs):
        """Transform raw input data into forcing variables for one timestep.

        Parameters
        ----------
        date : datetime-like, optional
            Target date/time used when loading the input dataset if ``ds`` is
            not provided.
        ds : xarray.Dataset, optional
            Already opened input dataset. If ``None`` it is loaded via
            :meth:`get_input_dataset`.
        outpath : str, optional
            Directory where output should be written. If ``None`` the
            preprocessed dataset is returned instead of writing.
        initial : bool or None, optional
            Flag passed through to remapping indicating initial-condition
            specific processing. If ``None`` subclasses may infer for first
            timestep.
        **kwargs
            Additional keyword arguments forwarded to the ``remap`` function.

        Returns
        -------
        xarray.Dataset or str
            In-memory forcing dataset (if ``outpath`` is ``None``) or path to
            the written file.
        """
        if ds is None:
            ds = self.get_input_dataset(date=date, initial=initial)
            if ds is None:
                warn(f"No input dataset available for date {date}")
        ads = self.remap(
            ds,
            domain_info=self.domain_info,
            vc=self.vc,
            surflib=self.surflib,
            initial=initial,
            **kwargs,
        )
        if outpath is None:
            return ads
        return self.write(ads, outpath)

    def run(
        self,
        start,
        end=None,
        freq="6h",
        outpath=None,
        compute=False,
        write=True,
        initial=None,
        **kwargs,
    ):
        """Batch preprocess a sequence of timesteps.

        Parameters
        ----------
        start : str or datetime-like
            Start date/time.
        end : str or datetime-like, optional
            Inclusive end date/time. If ``None`` only ``start`` is processed.
        freq : str, optional
            Timestep frequency passed to :func:`datelist` (default ``"6h"``).
        outpath : str, optional
            Output path template. Falls back to ``self.outpath`` if ``None``.
        compute : bool, optional
            If ``True`` triggers immediate Dask computation; otherwise delayed
            objects are returned.
        write : bool, optional
            If ``True`` datasets are written to disk; otherwise returned in
            memory.
        initial : bool or None, optional
            Explicit initial-condition flag. If ``None`` the first timestep is
            marked initial.
        **kwargs
            Extra keyword arguments forwarded to :meth:`preprocess` / ``remap``.

        Returns
        -------
        list
            List of ``xarray.Dataset`` objects or file paths depending on
            ``write`` / ``compute`` flags.
        """
        if outpath is None:
            outpath = self.outpath
        if end is None:
            end = start
        dates = datelist(start, end, freq=freq)
        if outpath is not None:
            outpath = outpath.format(date=dates[0])
            ensure_dir(outpath)
        print(f"writing to {outpath}")
        result = [
            self.preprocess(
                date=date,
                outpath=outpath if write is True else None,
                initial=initial or (initial is None and i == 0),
                **kwargs,
            )
            for i, date in enumerate(dates)
        ]
        if compute:
            result = dask.compute(*result)
        if write is True and compute is True:
            self._clean()
        return result


class CFPreprocessor(Preprocessor):
    """Preprocessor for CF-compliant GCM input datasets.

    Extends :class:`Preprocessor` by constructing a CF-style multi-variable
    input accessor (``gfile``) used to retrieve the raw dataset for each
    timestep.

    Parameters
    ----------
    expid : str
        Experiment identifier.
    surflib : str
        Path to surface library file.
    domain, vc, outpath, scratch : See :class:`Preprocessor`.
    input_data : dict, optional
        Specification of variables / paths needed by :func:`get_gcm_gfile`.
    **kwargs
        Additional keyword arguments passed through (reserved for future use).

    Attributes
    ----------
    gfile : object
        Accessor providing ``gfile(date)`` to obtain the raw input dataset.
    input_data : dict or None
        Original input specification.
    """

    def __init__(
        self,
        expid,
        surflib,
        domain=None,
        vc="vc_49lev",
        outpath=None,
        scratch=None,
        input_data=None,
        **kwargs,
    ):
        super().__init__(
            expid, surflib, domain=domain, vc=vc, outpath=outpath, scratch=scratch
        )
        self.input_data = input_data
        self.gfile = get_gcm_gfile(self.scratch.name, **self.input_data)
        self.remap = remap

    def get_input_dataset(self, date, initial=False, **kwargs):
        """
        Get the input dataset for a given date.

        Parameters
        ----------
        date : datetime-like
            Date for the input dataset.

        Returns
        -------
        xarray.Dataset
            Input dataset.
        """
        return self.gfile.gfile(date, **kwargs).load()


class RemoPreprocessor(Preprocessor):
    """Preprocessor for REMO model output used as nesting input.

    Provides logic to locate existing REMO NetCDF files and remap them for a
    subsequent nested run.

    Parameters
    ----------
    expid : str
        Target experiment identifier for new forcing.
    surflib : str
        Path to surface library file.
    domain, vc, outpath, scratch : See :class:`Preprocessor`.
    input_data : dict, optional
        Dictionary containing at least ``{"path": <input_dir>, "expid": <source_expid>}``.

    Attributes
    ----------
    inpath : str
        Directory containing source REMO files.
    inexp : str
        Source experiment id inside filenames.
    filename_pattern : str
        Python format string used to construct input filenames.
    """

    def __init__(
        self,
        expid,
        surflib,
        domain=None,
        vc="vc_49lev_nh_pt2000",
        outpath=None,
        scratch=None,
        input_data=None,
    ):
        super().__init__(
            expid, surflib, domain=domain, vc=vc, outpath=outpath, scratch=scratch
        )
        self.input_data = input_data
        self.remap = remap_remo
        self.inpath = input_data["path"]
        self.inexp = input_data["expid"]
        self.filename_pattern = "e{expid}t{date:%Y%m%d%H}.nc"

    def get_filename(self, date):
        """Construct the path of the source REMO file for a given date.

        Parameters
        ----------
        date : datetime-like
            Timestep whose file should be referenced.

        Returns
        -------
        str
            Absolute path to the expected REMO input file.
        """
        return os.path.join(
            self.inpath,
            self.filename_pattern.format(expid=self.inexp, date=date),
        )

    def open_remo_dataset(self, filename):
        """
        Open a REMO dataset and parse its dates.

        Parameters
        ----------
        filename : str
            Path to the REMO dataset file.

        Returns
        -------
        xarray.Dataset
            Parsed REMO dataset.
        """
        return parse_dates(xr.open_dataset(filename), use_cftime=True)

    def get_input_dataset(self, date, initial=False):
        """
        Get the input dataset for a given date.

        Parameters
        ----------
        date : datetime-like
            Date for the input dataset.

        Returns
        -------
        xarray.Dataset
            Input dataset.
        """
        filename = self.get_filename(date)
        if not os.path.isfile(filename):
            warn(f"Expected input file not found: {filename}")
            return None
        return self.open_remo_dataset(filename)


class ERA5Preprocessor(Preprocessor):
    """Preprocessor for ERA5 reanalysis data.

    Downloads / constructs hourly ERA5 forcing files locally via
    :func:`era5_gfile_from_dkrz` and remaps them onto the target domain.

    Parameters
    ----------
    expid : str
        Experiment identifier.
    surflib : str
        Path to surface library file.
    domain, vc, outpath, scratch : See :class:`Preprocessor`.
    input_data : dict, optional
        Reserved for future configuration hooks (currently unused).

    Attributes
    ----------
    input_data : dict or None
        Configuration passed by caller.
    """

    def __init__(
        self,
        expid,
        surflib,
        domain=None,
        vc="vc_49lev",
        outpath=None,
        scratch=None,
        input_data=None,
    ):
        super().__init__(
            expid, surflib, domain=domain, vc=vc, outpath=outpath, scratch=scratch
        )
        self.input_data = input_data

    def get_input_dataset(self, date=None, filename=None, initial=False):
        """
        Get the input dataset for a given date.

        Parameters
        ----------
        date : datetime-like
            Date for the input dataset.

        Returns
        -------
        xarray.Dataset
            Input dataset.
        """
        if filename is None:
            filename = era5_gfile_from_dkrz(date, self.scratch.name, add_soil=initial)
            print(f"created: {filename}")
        # logger.debug(f"created: {filename}")
        ds = xr.open_dataset(filename).load()
        return get_gcm_dataset(ds)


class CloudPreprocessor(Preprocessor):
    """Preprocessor sourcing CMIP6-like input from public cloud catalogs.

    Uses an intake-esm catalog to assemble required atmospheric and oceanic
    variables stored in cloud object storage and merges them into a single
    dataset used for forcing generation.

    Parameters
    ----------
    expid : str
        Experiment identifier.
    surflib : str
        Path to surface library file.
    domain, vc, outpath, scratch : See :class:`Preprocessor`.
    input_data : dict, optional
        Search parameters overriding catalog defaults (e.g. ``{"source_id": "MPI-ESM1-2-HR"}``).
    url : str, optional
        Shortcut key (``"gc"`` / ``"aws"``) or full URL to pangeo CMIP6 JSON catalog.

    Attributes
    ----------
    cat : intake.catalog
        Filtered intake-esm catalog.
    dsets : dict[str, xarray.Dataset]
        Individual variable datasets prior to merge.
    gcm : xarray.Dataset
        Merged multi-variable atmospheric dataset.
    tos : xarray.Dataset
        Sea surface temperature dataset.
    gfile : xarray.Dataset
        Convenience merged dataset providing ``sel(time=...)`` access.
    """

    atm = ["ta", "ua", "va", "hus", "sftlf", "orog"]
    ocn = ["tos"]

    def __init__(
        self,
        expid,
        surflib,
        domain=None,
        vc="vc_49lev",
        outpath=None,
        scratch=None,
        input_data=None,
        url=None,
    ):
        super().__init__(
            expid, surflib, domain=domain, vc=vc, outpath=outpath, scratch=scratch
        )
        self.input_data = input_data
        if url is None:
            url = "gc"
        if url in cloud_urls.keys():
            url = cloud_urls[url]
        self.url = url

        intake_default_kwargs = {
            "experiment_id": "historical",
            "member_id": "r1i1p1f1",
            "activity_id": "CMIP",
            #  "table_id": ["6hrLev", "fx", "Oday"],
            # "variable_id": self.atm + self.ocn,
        }
        intake_kwargs = intake_default_kwargs | input_data
        print(f"opening: {self.url}")
        self.cat = get_pangeo_catalog(self.url, **intake_kwargs)
        print("opening datasets...")
        self.dsets = get_dsets(self.cat.df[self.cat.df.variable_id.isin(self.atm)])
        print("merging datasets...")

        gcm = xr.merge(merge_dsets(self.dsets), join="override")

        tos = open_zstore(self.cat.df[self.cat.df.variable_id == "tos"].iloc[0].zstore)
        tos = tos.convert_calendar(tos.time.dt.calendar, use_cftime=True)
        self.tos = tos
        self.gcm = gcm
        self.gfile = get_gcm_dataset(gcm, tos=tos)

    def get_input_dataset(self, date, initial=False):
        """
        Get the input dataset for a given date.

        Parameters
        ----------
        date : datetime-like
            Date for the input dataset.

        Returns
        -------
        xarray.Dataset
            Input dataset.
        """
        return self.gfile.sel(time=date).compute()
