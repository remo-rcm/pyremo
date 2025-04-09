import os
import tempfile

import dask
import xarray as xr

import pyremo as pr
from ..remo_ds import update_meta_info, parse_dates
from .era5 import era5_gfile_from_dkrz
from .utils import datelist, ensure_dir, write_forcing_file
from .remapping import remap, remap_remo
from .cf import get_gcm_dataset, get_gcm_gfile


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
    """
    Preprocessor class for preparing input data.

    Parameters
    ----------
    expid : str
        Experiment ID.
    domain : dict or str
        Domain information.
    vc : str or dict
        Vertical coordinate information.
    surflib : str
        Path to the surface library.
    outpath : str, optional
        Output path, by default None.
    scratch : str, optional
        Scratch directory, by default None.
    """

    def __init__(
        self, expid, surflib, domain=None, vc="vc_49lev", outpath=None, scratch=None
    ):
        if scratch is None:
            scratch = os.environ["SCRATCH"]
        self.scratch = tempfile.TemporaryDirectory(dir=scratch)
        # logger.debug(f"scratch: {self.scratch.name}")
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
        """
        Clean up the scratch directory.
        """
        # logger.debug(f"cleaning up: {self.scratch.name}")
        self.scratch.cleanup()

    def _init_domain_info(self, domain=None):
        if isinstance(domain, dict):
            self.domain_info = domain
        elif isinstance(domain, str):
            self.domain_info = pr.domain_info(domain)
        else:
            domain_info = self.surflib.cx.info()
            domain_info["nlon"] -= 2
            domain_info["nlat"] -= 2
            domain_info["ll_lon"] += domain_info["dlon"]
            domain_info["ll_lat"] += domain_info["dlat"]
            self.domain_info = domain_info

    def write(self, ds, outpath):
        """
        Write the dataset to a forcing file.

        Parameters
        ----------
        ds : xarray.Dataset
            Dataset to write.
        outpath : str
            Output path.

        Returns
        -------
        str
            Path to the written file.
        """
        return write_forcing_file(ds, path=outpath, expid=self.expid)

    @dask.delayed
    def preprocess(self, date=None, ds=None, outpath=None, initial=False, **kwargs):
        """
        Preprocess the dataset for a given date.

        Parameters
        ----------
        date : datetime-like
            Date for preprocessing.
        outpath : str, optional
            Output path, by default None.

        Returns
        -------
        xarray.Dataset or str
            Preprocessed dataset or path to the written file.
        """
        if ds is None:
            ds = self.get_input_dataset(date=date)
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
        """
        Run the preprocessing for a given date range.

        Parameters
        ----------
        start : str or datetime-like
            Start date.
        end : str or datetime-like, optional
            End date, by default None.
        freq : str, optional
            Frequency, by default "6h".
        outpath : str, optional
            Output path, by default None.
        compute : bool, optional
            Whether to compute the results, by default False.
        write : bool, optional
            Whether to write the results, by default True.

        Returns
        -------
        list
            List of preprocessed datasets or paths to the written files.
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
    """
    CFPreprocessor class for preparing input data from CF-compliant datasets.

    Parameters
    ----------
    input_data : dict
        Input data information.
    *args : tuple
        Additional arguments.
    **kwargs : dict
        Additional keyword arguments.
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

    def get_input_dataset(self, date):
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
        return self.gfile.gfile(date).load()


class RemoPreprocessor(Preprocessor):
    """
    RemoPreprocessor class for preparing double nesting input data.

    Parameters
    ----------
    expid : str
        Experiment ID.
    surflib : str
        Path to the surface library.
    domain : dict or str, optional
        Domain information, by default None.
    vc : str or dict, optional
        Vertical coordinate information, by default "vc_49lev_nh_pt2000".
    outpath : str, optional
        Output path, by default None.
    scratch : str, optional
        Scratch directory, by default None.
    input_data : dict, optional
        Input data information, by default None.
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

    def get_input_dataset(self, date):
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
        return self.open_remo_dataset(filename)


class ERA5Preprocessor(Preprocessor):
    """
    ERA5Preprocessor class for preparing input data from ERA5 datasets.

    Parameters
    ----------
    input_data : dict
        Input data information.
    *args : tuple
        Additional arguments.
    **kwargs : dict
        Additional keyword arguments.
    """

    def __init__(
        self, expid, surflib, domain=None, vc="vc_49lev", outpath=None, scratch=None
    ):
        super().__init__(
            expid, surflib, domain=domain, vc=vc, outpath=outpath, scratch=scratch
        )

    def get_input_dataset(self, date=None, filename=None):
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
            filename = era5_gfile_from_dkrz(date, self.scratch.name)
            print(f"created: {filename}")
        # logger.debug(f"created: {filename}")
        ds = xr.open_dataset(filename).load()
        return get_gcm_dataset(ds)


class CloudPreprocessor(Preprocessor):

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

    def get_input_dataset(self, date):
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
