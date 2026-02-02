"""Cmorizer

ERA5 cmorizer is mostly based on the official ECMWF documentation of converting GRIB to NetCDF:

https://confluence.ecmwf.int/display/OIFS/How+to+convert+GRIB+to+netCDF

"""

import os
import subprocess
from os import path as op
from pprint import pprint
from warnings import warn
from importlib.resources import files
import copy

import pandas as pd
import xarray as xr
from cdo import Cdo

from ..utils import read_yaml

# path and file templates at DKRZ
# see https://docs.dkrz.de/doc/dataservices/finding_and_accessing_data/era_data/#file-and-directory-names

dkrz_template = {
    "path_template": "/pool/data/ERA5/{era_id}/{level_type}/{dataType}/{frequency}/{code:03d}",
    "file_template": "{era_id}{level_type}{typeid}_{frequency}_{date}_{code:03d}.grb",
}

era5_params_file = files("pyremo.preproc").joinpath("era5-params.yaml")
era5_grid_file = files("pyremo.preproc").joinpath("grid.txt")
ecmwf_params_file = files("pyremo.preproc").joinpath("ecmwf_128.csv")

era5_params_config = read_yaml(era5_params_file)

era5_params = {
    k: era5_params_config.get("defaults", {}) | v
    for k, v in era5_params_config["parameters"].items()
}

ecmwf_params = pd.read_csv(ecmwf_params_file)


def get_output_filename(date, expid, path=None):
    """Build standard REMO gfile output filename.

    Parameters
    ----------
    date : str or datetime-like
        Datetime used in the filename; formatted as ``YYYYMMDDHH``.
    expid : str
        Experiment identifier to include in the filename (e.g., ``"000000"``).
    path : str, optional
        Directory where the file will be written. Defaults to current directory.

    Returns
    -------
    str
        Full path to the output file in the form ``g{expid}a{YYYYMMDDHH}.nc``.
    """
    if path is None:
        path = "./"
    date = pd.to_datetime(date).strftime("%Y%m%d%H")
    return op.join(path, f"g{expid}a{date}.nc")


def params_by_code(params, code):
    """Return parameter configuration matching a GRIB code.

    Parameters
    ----------
    params : dict
        Mapping of variable name to parameter configuration dictionaries.
    code : int
        GRIB parameter code to search for.

    Returns
    -------
    dict or None
        The matching parameter configuration merged with ``{"rename": varname}``,
        or ``None`` if no match is found.
    """
    for k, v in params.items():
        c = v.get("code", -999)
        if code == c:
            return v | {"rename": k}
    return None


def renames_by_code(ds, params):
    """Build a rename mapping based on GRIB codes in a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset whose variables include a ``code`` attribute.
    params : dict
        Parameter configuration used to map GRIB codes to target names.

    Returns
    -------
    dict
        Mapping from original variable names to target names.
    """
    renames = {}
    for var in ds.data_vars:
        code = ds[var].attrs.get("code", -999)
        param = params_by_code(params, code)
        if param:
            renames[var] = param["rename"]
    return renames


def get_params(config):
    """Expand parameter configuration with defaults.

    Parameters
    ----------
    config : dict
        Configuration containing a ``defaults`` dict and a ``parameters`` mapping.

    Returns
    -------
    dict
        Mapping of parameter keys to dictionaries with defaults applied.
    """
    defaults = config.get("defaults", {})
    return {k: defaults | v for k, v in config["parameters"].items()}


def check_search(result):
    """Validate an intake search result and return a single path.

    Parameters
    ----------
    result : pandas.DataFrame
        DataFrame returned by ``catalog.search(...).df``.

    Returns
    -------
    str or None
        The selected path if a unique result exists, otherwise ``None``.
    """
    if len(result) == 0:
        warn(f"nothing found: {len(result)}")
        return None
    if len(result) > 1:
        pass
        # warn(f"result not unique: {len(result)}, {result}")
    return result.iloc[0].path


def get_file_from_intake(cat, params, date=None):
    """Query an intake catalog for a single parameter entry.

    Parameters
    ----------
    cat : intake.esm.EsmDatastore
        Intake catalog to query.
    params : dict
        Search parameters describing the ERA5 entry (e.g., code, level_type).
    date : str or datetime-like, optional
        Date used to filter entries; coerced to ``YYYY-MM-DD``. If the
        parameter frequency is invariant, ``INVARIANT`` is used.

    Returns
    -------
    intake.catalog.local.LocalCatalogEntry
        The search result object; call ``.df`` to access the underlying table.
    """
    if date is not None:
        date = pd.to_datetime(date).strftime("%Y-%m-%d")
    freq = params.get("frequency", None)
    if freq is None or freq == "invariant":
        date = "INVARIANT"
    return cat.search(**params, validation_date=date)


def get_files_from_intake(cat, params, date=None):
    """Collect file paths from an intake catalog for multiple variables.

    Parameters
    ----------
    cat : intake.esm.EsmDatastore
        Intake catalog to query.
    params : dict
        Mapping from variable name to search parameters.
    date : str or datetime-like, optional
        Date used to filter entries.

    Returns
    -------
    dict
        Mapping from variable name to resolved file path (or ``None`` if not found).
    """
    files = {}
    for k, v in params.items():
        f = check_search(get_file_from_intake(cat, params=v, date=date).df)
        if not f:
            warn(f"no result for {k} --> params: {v}")
        files[k] = f
    return files


def get_file_from_template(
    date,
    era_id,
    frequency,
    dataType,
    code,
    level_type,
    template=None,
    **kwargs,
):
    """Derive a filename from the DKRZ template.

    Parameters
    ----------
    date : str or datetime-like
        Date to include in the path/filename; becomes ``YYYY-MM-DD`` unless invariant.
    era_id : str
        ERA dataset identifier (e.g., ``"E5"`` or ``"E1"``).
    frequency : {"hourly", "daily", "monthly", "invariant"}
        Data frequency.
    dataType : {"an", "fc"}
        Data type used to derive the type id in the filename.
    code : int
        GRIB parameter code.
    level_type : {"model_level", "surface"}
        Level type selector.
    template : dict, optional
        Template dictionary with ``path_template`` and ``file_template``. Defaults to
        the built-in DKRZ template.

    Returns
    -------
    str
        Absolute template-expanded path to the GRIB file.
    """
    if template is None:
        template = dkrz_template

    path_template = template["path_template"]
    file_template = template["file_template"]

    lt = {
        "model_level": "ml",
        "surface": "sf",
    }
    freqs = {
        "hourly": "1H",
        "daily": "1D",
        "monthly": "1M",
        "invariant": "IV",
    }
    typeids = {
        "an": "00",
        "fc": "12",
    }

    level_type = lt.get(level_type, level_type)
    frequency = freqs.get(frequency, frequency)
    typeid = typeids.get(dataType, dataType)
    if frequency == "IV":
        date = "INVARIANT"
    else:
        date = pd.to_datetime(date).strftime("%Y-%m-%d")

    filepath_template = op.join(path_template, file_template)

    return filepath_template.format(
        date=date,
        era_id=era_id,
        frequency=frequency,
        dataType=dataType,
        code=code,
        level_type=level_type,
        typeid=typeid,
    )


def get_files_from_template(params, date, template=None):
    """Resolve GRIB file paths for multiple variables from a template.

    Parameters
    ----------
    params : dict
        Mapping from variable name to template fields (see
        ``get_file_from_template``).
    date : str or datetime-like
        Desired date.
    template : dict, optional
        Template dictionary with ``path_template`` and ``file_template``.

    Returns
    -------
    dict
        Mapping from variable name to file path.
    """
    files = {}
    for k, v in params.items():
        f = get_file_from_template(date=date, template=template, **v)
        if not op.isfile(f) and v.get("era_id") == "E1":
            warn(f"file not found: {f}, will try E5")
            v["era_id"] = "E5"
            f = get_file_from_template(date=date, template=template, **v)
        files[k] = f
    return files


class ERA5:
    """
    Class for cmorizing original ERA5 GRIB data.

    Notes
    -----
    The cmorizer class mostly works with the intake catalog provided by DKRZ.

    References
    ----------
    Please refer to the DKRZ data pool documentation: https://docs.dkrz.de/doc/dataservices/finding_and_accessing_data/era_data

    """

    dynamic = ["ta", "hus", "ps", "tos", "sic", "clw", "snd"]
    wind = ["svo", "sd"]
    fx = ["orog", "sftlf"]
    soil_vars = [
        "tsl1",
        "tsl2",
        "tsl3",
        "tsl4",
        "swvl1",
        "swvl2",
        "swvl3",
        "swvl4",
        "tsn",
        "src",
        "skt",
        "tos",
        "sic",
        "skt",
        "snd",
    ]

    chunks = {}
    options = "-f nc4"

    def __init__(
        self, params=None, cat=None, gridfile=None, scratch=None, template=None
    ):
        """Create a new ERA5 cmorizer instance.

        Parameters
        ----------
        params : dict, optional
            Parameter configuration mapping (defaults to built-in ERA5 parameters).
        cat : intake.esm.EsmDatastore or str, optional
            Intake catalog or path/URL to an ESM datastore.
        gridfile : str, optional
            Path to CDO grid description file. Defaults to package resource.
        scratch : str, optional
            Temporary directory for CDO operations. Defaults to ``$SCRATCH`` or ``."``.
        template : dict, optional
            File path template mapping. Defaults to the DKRZ template.
        """
        if isinstance(cat, str):
            import intake

            self.cat = intake.open_esm_datastore(cat)
        else:
            self.cat = cat
        if scratch is None:
            scratch = os.environ.get("SCRATCH", "./")
        self.scratch = scratch
        if params is None:
            params = era5_params
        self.params = params
        if gridfile is None:
            gridfile = era5_grid_file
        self.gridfile = gridfile
        if template is None:
            template = dkrz_template
        self.template = template
        self.cdo = Cdo(tempdir=scratch)

    def _get_files(self, date, variables=None):
        """Resolve input files for a set of variables and a given date.

        Parameters
        ----------
        date : str or datetime-like
            Target date for selection.
        variables : list of str, optional
            Variable names to include. Defaults to dynamic, fx, and wind variables.

        Returns
        -------
        dict
            Mapping from variable name to file path.
        """
        if variables is None:
            variables = self.dynamic + self.fx + self.wind
        if self.cat:
            return get_files_from_intake(
                self.cat,
                {k: v for k, v in self.params.items() if k in variables},
                date,
            )
        else:
            return get_files_from_template(
                date=date,
                template=self.template,
                params={k: v for k, v in self.params.items() if k in variables},
            )

    def _seldate(self, filename, date):
        """Build a CDO command to select a single date from a file."""
        return f"-seldate,{date} {filename}"
        # return self.cdo.seldate(date, input=filename)

    def _seldates(self, filenames, date):
        """Apply ``_seldate`` to all non-invariant files in a mapping."""
        return {
            v: (
                self._seldate(f, date)
                if self.params[v].get("frequency") != "invariant"
                else f
            )
            for v, f in filenames.items()
        }

    def _to_regulars(self, filenames, gridtypes):
        """Convert multiple inputs to regular Gaussian grid commands."""
        return {
            v: self._to_regular(f, gridtype=gridtypes[v], setname=v)
            for v, f in filenames.items()
            if v in self.dynamic + self.fx + self.soil_vars
        }

    def _get_gridtypes(self, filenames):
        """Determine grid types for all input files using CDO griddes."""
        return {v: self._gridtype(f) for v, f in filenames.items()}

    def _open_dsets(self, filenames):
        """Open a mapping of files as chunked ``xarray.Dataset`` objects."""
        dsets = {}
        for v, f in filenames.items():
            ds = xr.open_dataset(f, chunks=self.chunks)
            if v in self.fx:
                # squeeze out time
                ds = ds.squeeze(drop=True)
            dsets[v] = ds
        return dsets

    def _griddes(self, filename):
        """Return parsed CDO griddes information as a dictionary."""
        griddes = self.cdo.griddes(input=filename)
        return {
            entry.split("=")[0].strip(): entry.split("=")[1].strip()
            for entry in griddes
            if "=" in entry
        }

    def _gridtype(self, filename):
        """Extract grid type from CDO griddes output."""
        return self._griddes(filename)["gridtype"]

    def _to_regular(self, filename, gridtype=None, setname="", table="ecmwf"):
        """Convert ECMWF GRIB to regular Gaussian NetCDF (CDO command).

        CDO is used depending on the detected ``gridtype``:

        - ``gaussian_reduced`` → ``setgridtype,regular``
        - ``spectral`` → ``sp2gpl``
        - ``gaussian`` → passthrough

        This follows the ECMWF ERA5 recommendation.
        """
        if table is None:
            table = ""
        # from cdo import Cdo
        # cdo = Cdo(tempdir=scratch)
        if gridtype is None:
            gridtype = self._gridtype(filename)
        # options = f"-f nc4 -t {table}"
        if setname:
            setname = f"-setname,{setname}"  # {filename}"
        if gridtype == "gaussian_reduced":
            gaussian = "-setgridtype,regular"
        elif gridtype == "spectral":
            gaussian = "-sp2gpl"
        elif gridtype == "gaussian":
            gaussian = ""
        else:
            raise Exception(
                "unknown grid type for conversion to regular grid: {}".format(gridtype)
            )
        command = f"{setname} {gaussian} {filename}"
        return command

    def _compute_wind(self, vort, div):
        """Compute wind from vorticity and divergence (CDO command)."""
        return f"-chname,u,ua,v,va -dv2uvl -merge [ {vort} {div} ]"

    def gfile(self, date, path=None, expid=None, filename=None, add_soil=False):
        """Create an ERA5 gfile NetCDF containing required variables.

        Converts ERA5 GRIB inputs to a regular Gaussian grid and merges all
        required variables into a single file for REMO preprocessing.

        Parameters
        ----------
        date : str or datetime-like
            Date (ISO 8601) for which variables should be converted.
        path : str, optional
            Output directory for the gfile.
        expid : str, optional
            Experiment id for filename templating. Defaults to ``"000000"``.
        filename : str, optional
            Full output filename. If omitted, it is composed from ``path`` and
            ``expid``.
        add_soil : bool, optional
            If True, include additional soil variables.

        Returns
        -------
        str
            The output filename.
        """
        if expid is None:
            expid = "000000"
        if filename is None:
            filename = get_output_filename(date, expid, path)

        variables = self.dynamic + self.fx + self.wind
        if add_soil is True:
            variables += self.soil_vars

        print(f"output filename: {filename}")
        # gridfile = "/work/ch0636/g300046/remo/era5-cmor/notebooks/grid.txt"

        print("getting files...")
        files = self._get_files(date, variables=variables)
        pprint(f"using files: \n{files}")
        print("getting gridtypes...")
        gridtypes = self._get_gridtypes(files)
        print("selecting dates...")
        seldates = self._seldates(files, date)
        print("convert to regular grid...")
        regulars = self._to_regulars(seldates, gridtypes)
        print("computing wind...")
        wind = self._compute_wind(seldates["svo"], seldates["sd"])

        merge = (
            f"-setgrid,{self.gridfile} -merge [ "
            + " ".join(list(regulars.values()) + [wind])
            + " ]"
        )
        call = f"cdo {self.options} invertlev -invertlat {merge} {filename}"
        print(f"execute: {call}")

        subprocess.run(
            call.split(),
            check=True,
            shell=False,
        )
        # stdout, stderr = process.communicate()

        return filename

    def get_soil(self, date, filename=None):
        """Create a NetCDF with ERA5 soil-related variables.

        Parameters
        ----------
        date : str or datetime-like
            Target date.
        filename : str, optional
            Output filename. Defaults to ``$SCRATCH/era5_soil_data.nc``.

        Returns
        -------
        str
            The output filename.
        """
        if not isinstance(date, str):
            date = date.strftime("%Y-%m-%dT%H:%M:%S")
        if filename is None:
            filename = op.join(self.scratch, "era5_soil_data.nc")
        files = self._get_files(date, self.soil_vars + self.fx)
        gridtypes = self._get_gridtypes(files)
        seldates = self._seldates(files, date)
        regulars = self._to_regulars(seldates, gridtypes)
        merge = (
            f"-setgrid,{self.gridfile} -merge [ "
            + " ".join(list(regulars.values()))
            + " ]"
        )
        call = f"cdo {self.options} invertlat {merge} {filename}"
        print(f"execute: {call}")
        subprocess.run(
            call.split(),
            check=True,
            shell=False,
        )

        return filename


def era5_from_gcloud(time=None, chunks="auto", hfreq=6):
    """
    Create ERA5 input for preprocessing from Google Cloud.

    Parameters
    ----------
    time : str or datetime-like, optional
        The time range to select. Defaults to None.
    chunks : str or dict, optional
        The chunk size for Dask. Defaults to "auto".
    hfreq : int, optional
        The hourly frequency to select. Defaults to 6.

    Returns
    -------
    xarray.Dataset
        The selected ERA5 dataset.
    """

    # era5_to_cmip_cf = {
    #     "temperature": "ta",
    #     "surface_pressure": "ps",
    #     "u_component_of_wind": "ua",
    #     "v_component_of_wind": "va",
    #     "specific_humidity": "hus",
    #     "specific_cloud_liquid_water_content": "clw",
    #     "geopotential_at_surface": "orog",
    #     "land_sea_mask": "sftlf",
    #     "sea_ice_cover": "sic",
    #     "snow_depth": "snd",
    #     "sea_surface_temperature": "tos",
    # }

    era5_to_cmip_cf = {v["gc_name"]: k for k, v in era5_params.items()}

    if chunks == "auto":
        chunks = {"time": 48}

    ds = xr.open_zarr(
        "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
        chunks=chunks,
        storage_options=dict(token="anon"),
        consolidated=True,
    )

    ar_full_37_1h = ds.sel(
        time=slice(ds.attrs["valid_time_start"], ds.attrs["valid_time_stop"], hfreq)
    )

    ds = xr.open_zarr(
        "gs://gcp-public-data-arco-era5/ar/model-level-1h-0p25deg.zarr-v1",
        chunks=chunks,
        storage_options=dict(token="anon"),
        consolidated=True,
    )
    ar_native_vertical_grid_data = ds.sel(
        time=slice(ds.attrs["valid_time_start"], ds.attrs["valid_time_stop"], hfreq)
    )

    if time:
        ar_full_37_1h = ar_full_37_1h.sel(time=time)
        ar_native_vertical_grid_data = ar_native_vertical_grid_data.sel(time=time)

    level_vars = [
        v for v in era5_to_cmip_cf.keys() if v in ar_native_vertical_grid_data.data_vars
    ]
    surface_vars = [v for v in era5_to_cmip_cf.keys() if v in ar_full_37_1h.data_vars]

    ds = xr.merge(
        [
            ar_native_vertical_grid_data[level_vars],
            ar_full_37_1h[surface_vars].drop_dims("level"),
        ]
    )

    ds = ds.rename(era5_to_cmip_cf)

    # if freq:
    #    ds = ds.isel(time=slice(0, None, freq))

    for v in ["orog", "sftlf"]:
        ds[v] = ds[v].isel(time=0).squeeze(drop=True)

    pv = ds.ta.GRIB_pv
    n = len(pv)
    hyai, hybi = pv[0 : n // 2], pv[n // 2 : n]

    ds["akgm"] = xr.DataArray(hyai, dims="nhyi")
    ds["bkgm"] = xr.DataArray(hybi, dims="nhyi")

    ds = ds.rename(hybrid="lev", latitude="lat", longitude="lon").cf.guess_coord_axis(
        verbose=True
    )
    ds.lev.attrs["positive"] = "down"
    # invert lat
    ds = ds.reindex(lat=ds.lat[::-1])

    return ds


def era5_gfile_from_dkrz(date, path, expid="000000", add_soil=False):
    """
    Compute the gfile for a given date.

    Parameters
    ----------
    date : datetime-like
        The date for which to compute the gfile.
    path : str
        The path to save the gfile.
    expid : str
        The experiment ID.

    Returns
    -------
    str
        The path to the computed gfile.
    """
    params = copy.deepcopy(era5_params)

    if date.year in range(2000, 2007):
        for k, v in params.items():
            # if v["level_type"] == "model_level":
            v["era_id"] = "E1"

    era5 = ERA5(params, gridfile=era5_grid_file)
    return era5.gfile(
        date.strftime("%Y-%m-%dT%H:%M:%S"), path, expid, add_soil=add_soil
    )
