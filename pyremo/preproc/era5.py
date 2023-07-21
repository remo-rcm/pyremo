"""Cmorizer

ERA5 cmorizer is mostly based on the official ECMWF documentation of converting GRIB to NetCDF:

https://confluence.ecmwf.int/display/OIFS/How+to+convert+GRIB+to+netCDF

"""

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


class ERA5:
    dynamic = ["ta", "hus", "ps", "tos", "sic", "clw", "snd"]
    wind = ["svo", "sd"]
    fx = ["orog", "sftlf"]
    chunks = {}

    def __init__(self, cat, params, scratch=None, show_cdo=True):
        if isinstance(cat, str):
            import intake

            self.cat = intake.open_esm_datastore(cat)
        else:
            self.cat = cat
        self.scratch = scratch
        self.params = params
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
        return self.cdo.seldate(date, input=filename)

    def _seldates(self, filenames, date):
        return {
            v: (
                self._seldate(f, date)
                if self.params[v].get("frequency") != "invariant"
                else f
            )
            for v, f in filenames.items()
        }

    def _to_netcdf(self, filenames):
        return {
            v: self._to_regular(f, setname=v)
            for v, f in filenames.items()
            if v in self.dynamic + self.fx
        }

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

    def _to_regular(self, filename, table="ecmwf", setname=""):
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
        gridtype = self._gridtype(filename)
        options = f"-f nc4 -t {table}"
        if setname:
            filename = f"--setname,{setname} {filename}"
        print(filename)
        print(gridtype)
        if gridtype == "gaussian_reduced":
            gaussian = self.cdo.setgridtype("regular", options=options, input=filename)
        elif gridtype == "spectral":
            gaussian = self.cdo.sp2gpl(options=options, input=filename)
        elif gridtype == "gaussian":
            gaussian = self.cdo.copy(options=options, input=filename)
        else:
            raise Exception(
                "unknown grid type for conversion to regular grid: {}".format(gridtype)
            )
        return self.cdo.invertlat(input=gaussian)

    def _compute_wind(self, vort, div):
        """compute wind from vorticity and divergence"""
        # ds = xr.merge([vort_ds, div_ds])
        merge = self.cdo.merge(input=[vort, div])
        uv = self.cdo.dv2uvl(options="-f nc4", input=merge)
        uv = self.cdo.invertlat(input=uv)
        return xr.open_dataset(uv, chunks=self.chunks).rename({"u": "ua", "v": "va"})

    def gfile(
        self,
        date,
        write=True,
    ):
        """Create an ERA5 gfile dataset.

        Main function to convert ERA5 grib data to a regular gaussian Dataset
        containing all variables required for REMO preprocessing.

        Parameters
        ----------
        dates : date or list of dates in ISO 8601 format.
            A single date or list of dates for which the variables should
            be converted.
        cf_meta : bool
            Rename variables to CF standard names.
        clean_coords: bool
            Drop time coordinate from coordinate variables.
        add_fx: bool
            Add static fields, e.g., orography and land sea mask.
        fx_date: str
            Date used for adding fx variables. If ``fx_date=None``,
            the default date `"1979-01-01T00:00:00"` is used. Defaults
            to ``None``.

        Returns
        -------
        Dataset

        """
        print("getting files...")
        files = self._get_files(date)
        print("selecting dates...")
        seldates = self._seldates(files, date)
        print("to netcdf..")
        regulars = self._to_netcdf(seldates)
        print("computing wind...")
        wind = self._compute_wind(seldates["svo"], seldates["sd"])
        print("merging...")
        dsets = self._open_dsets(regulars)

        gds = xr.merge(
            list(dsets.values()) + [wind], join="override", compat="override"
        )

        return gds
