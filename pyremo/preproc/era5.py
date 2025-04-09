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
    if path is None:
        path = "./"
    date = pd.to_datetime(date).strftime("%Y%m%d%H")
    return op.join(path, f"g{expid}a{date}.nc")


def params_by_code(params, code):
    for k, v in params.items():
        c = v.get("code", -999)
        if code == c:
            return v | {"rename": k}
    return None


def renames_by_code(ds, params):
    renames = {}
    for var in ds.data_vars:
        code = ds[var].attrs.get("code", -999)
        param = params_by_code(params, code)
        if param:
            renames[var] = param["rename"]
    return renames


def get_params(config):
    defaults = config.get("defaults", {})
    return {k: defaults | v for k, v in config["parameters"].items()}


def check_search(result):
    if len(result) == 0:
        warn(f"nothing found: {len(result)}")
        return None
    if len(result) > 1:
        pass
        # warn(f"result not unique: {len(result)}, {result}")
    return result.iloc[0].path


def get_file_from_intake(cat, params, date=None):
    """get filename entry from intake catalog"""
    if date is not None:
        date = pd.to_datetime(date).strftime("%Y-%m-%d")
    freq = params.get("frequency", None)
    if freq is None or freq == "invariant":
        date = "INVARIANT"
    return cat.search(**params, validation_date=date)


def get_files_from_intake(cat, params, date=None):
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
    """Derive filename from filename template

    Derives filename according to https://docs.dkrz.de/doc/dataservices/finding_and_accessing_data/era_data/#file-and-directory-names

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
    return op.join(path_template, file_template).format(
        date=date,
        era_id=era_id,
        frequency=frequency,
        dataType=dataType,
        code=code,
        level_type=level_type,
        typeid=typeid,
    )


def get_files_from_template(params, date, template=None):
    files = {}
    for k, v in params.items():
        f = get_file_from_template(date=date, template=template, **v)
        if not f:
            warn(f"no result for {k} --> params: {v}")
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
        return f"-seldate,{date} {filename}"
        # return self.cdo.seldate(date, input=filename)

    def _seldates(self, filenames, date):
        return {
            v: (
                self._seldate(f, date)
                if self.params[v].get("frequency") != "invariant"
                else f
            )
            for v, f in filenames.items()
        }

    def _to_regulars(self, filenames, gridtypes):
        return {
            v: self._to_regular(f, gridtype=gridtypes[v], setname=v)
            for v, f in filenames.items()
            if v in self.dynamic + self.fx + self.soil_vars
        }

    def _get_gridtypes(self, filenames):
        return {v: self._gridtype(f) for v, f in filenames.items()}

    def _open_dsets(self, filenames):
        dsets = {}
        for v, f in filenames.items():
            ds = xr.open_dataset(f, chunks=self.chunks)
            if v in self.fx:
                # squeeze out time
                ds = ds.squeeze(drop=True)
            dsets[v] = ds
        return dsets

    def _griddes(self, filename):
        griddes = self.cdo.griddes(input=filename)
        return {
            entry.split("=")[0].strip(): entry.split("=")[1].strip()
            for entry in griddes
            if "=" in entry
        }

    def _gridtype(self, filename):
        return self._griddes(filename)["gridtype"]

    def _to_regular(self, filename, gridtype=None, setname="", table="ecmwf"):
        """converts ecmwf spectral grib data to regular gaussian netcdf.

        cdo is used to convert ecmwf grid data to netcdf depending on the gridtype:
        For 'gaussian_reduced': cdo setgridtype,regular
            'spectral'        : cdo sp2gpl

        This follows the recommendation from the ECMWF Era5 Documentation.
        We also invert the latitudes to stick with cmor standard.

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
        """compute wind from vorticity and divergence"""
        return f"-chname,u,ua,v,va -dv2uvl -merge [ {vort} {div} ]"

    def gfile(self, date, path=None, expid=None, filename=None, add_soil=False):
        """Create an ERA5 gfile dataset.

        Main function to convert ERA5 grib data to a regular gaussian Dataset
        containing all variables required for REMO preprocessing.

        Parameters
        ----------
        date : date in ISO 8601 format.
            Date for which the variables should be converted.
        output: str
            Name of output file.
        path: str
            Output path for the gfile.
        expid: str
            Experiment id for the filenaming template.
        filename: str
            Filename including path for the output filename. If not provided,
            the filename will be created automatically from path and expid.

        Returns
        -------
        Output filename.

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

    def get_soil(self, date, path=None, expid=None, filename=None):
        files = self._get_files(date, self.soil_vars + self.fx)
        gridtypes = self._get_gridtypes(files)
        seldates = self._seldates(files, date)
        regulars = self._to_regulars(seldates, gridtypes)
        merge = (
            f"-setgrid,{self.gridfile} -merge [ "
            + " ".join(list(regulars.values()))
            + " ]"
        )
        call = f"cdo {self.options} invertlev -invertlat {merge} {filename}"
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


def era5_gfile_from_dkrz(date, path, expid="000000"):
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
    global era5_params

    params = era5_params.copy()
    if date.year in range(2000, 2007):
        for k, v in params.items():
            v["era_id"] = "E1"
    era5 = ERA5(params, gridfile=era5_grid_file)
    return era5.gfile(date.strftime("%Y-%m-%dT%H:%M:%S"), path, expid)
