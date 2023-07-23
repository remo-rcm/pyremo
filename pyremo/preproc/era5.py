"""Cmorizer

ERA5 cmorizer is mostly based on the official ECMWF documentation of converting GRIB to NetCDF:

https://confluence.ecmwf.int/display/OIFS/How+to+convert+GRIB+to+netCDF

"""

from os import path as op
from warnings import warn

import pandas as pd
import xarray as xr
from cdo import Cdo

# from ..utils import read_yaml

xr.set_options(keep_attrs=True)


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


def get_file(cat, params, date=None):
    if date is not None:
        date = pd.to_datetime(date).strftime("%Y-%m-%d")
    freq = params.get("frequency", None)
    if freq is None or freq == "invariant":
        date = "INVARIANT"
    return cat.search(**params, validation_date=date)


def get_files(cat, params, date=None):
    files = {}
    for k, v in params.items():
        f = check_search(get_file(cat, params=v, date=date).df)
        if not f:
            warn(f"no result for {k} --> params: {v}")
        files[k] = f
    return files


def get_filename(date, expid, path=None):
    if path is None:
        path = "./"
    date = pd.to_datetime(date).strftime("%Y%m%d%H")
    return op.join(path, f"g{expid}a{date}.nc")


class ERA5:
    dynamic = ["ta", "hus", "ps", "tos", "sic", "clw", "snd"]
    wind = ["svo", "sd"]
    fx = ["orog", "sftlf"]
    chunks = {}
    options = "-f nc4"

    def __init__(self, cat, params, gridfile=None, scratch=None):
        if isinstance(cat, str):
            import intake

            self.cat = intake.open_esm_datastore(cat)
        else:
            self.cat = cat
        self.scratch = scratch
        self.params = params
        self.gridfile = gridfile
        self.cdo = Cdo(tempdir=scratch)

    def _get_files(self, date):
        return get_files(
            self.cat,
            {
                k: v
                for k, v in self.params.items()
                if k in self.dynamic + self.fx + self.wind
            },
            date,
        )

    def _seldate(self, filename, date):
        return f"--seldate,{date} {filename}"
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
            if v in self.dynamic + self.fx
        }

    def _get_gridtypes(self, filenames):
        return {v: self._gridtype(f) for v, f in filenames.items()}

    def _to_netcdf(self, filenames):
        pass

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
            setname = f"--setname,{setname}"  # {filename}"
        if gridtype == "gaussian_reduced":
            gaussian = "--setgridtype,regular"
        elif gridtype == "spectral":
            gaussian = "--sp2gpl"
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
        return f"--chname,u,ua,v,va --dv2uvl --merge {vort} {div}"

    def gfile(self, date, path=None, expid=None, filename=None):
        """Create an ERA5 gfile dataset.

        Main function to convert ERA5 grib data to a regular gaussian Dataset
        containing all variables required for REMO preprocessing.

        Parameters
        ----------
        date : date in ISO 8601 format.
            Date for which the variables should be converted.
        output: str
            Name of output file.

        Returns
        -------
        Output filename.

        """
        if expid is None:
            expid = "000000"
        if filename is None:
            filename = get_filename(date, expid, path)
        print(f"filename: {filename}")
        # gridfile = "/work/ch0636/g300046/remo/era5-cmor/notebooks/grid.txt"

        print("getting files...")
        files = self._get_files(date)
        print("getting gridtypes...")
        gridtypes = self._get_gridtypes(files)
        print("selecting dates...")
        seldates = self._seldates(files, date)
        print("convert to regular grid...")
        regulars = self._to_regulars(seldates, gridtypes)
        print("computing wind...")
        wind = self._compute_wind(seldates["svo"], seldates["sd"])

        merge = f"--setgrid,{self.gridfile} --merge " + " ".join(
            list(regulars.values()) + [wind]
        )

        print(f"execute: {merge}")
        return self.cdo.invertlat(options=self.options, input=merge, output=filename)
