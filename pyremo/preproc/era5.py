"""Cmorizer

ERA5 cmorizer is mostly based on the official ECMWF documentation of converting GRIB to NetCDF:

https://confluence.ecmwf.int/display/OIFS/How+to+convert+GRIB+to+netCDF

"""
import os
import subprocess
import tempfile

import pandas as pd
import xarray as xr

try:
    from tqdm import tqdm
except Exception:

    def tqdm(x):
        return x


# try:
#     from dask.diagnostics import ProgressBar
#     pbar = ProgressBar()
#     pbar.register()
# except:
#     pass

#

from cordex import ecmwf as ectable

xr.set_options(keep_attrs=True)

# cdo_exe = "/sw/rhel6-x64/cdo/cdo-1.9.6-gcc64/bin/cdo"
cdo_exe = "cdo"

# os.path.join(os.environ['SCRATCH'], '.cdo_tmp')

tempfiles = []


def init_tempdir(dir=None):
    # global tempdir
    if dir is None:
        try:
            tempdir = os.path.join(os.environ["SCRATCH"], ".stoppistopp")
        except Exception:
            tempdir = tempfile.mkdtemp()
    else:
        tempdir = dir
    if not os.path.isdir(tempdir):
        os.makedirs(tempdir)
    return tempdir


varmap = {
    130: "ta",
    134: "ps",
    131: "ua",
    132: "va",
    133: "hus",
    34: "tos",
    31: "sic",
    129: "orog",
    172: "sftlf",
    139: "tsl1",
    170: "tsl2",
    183: "tsl3",
    238: "tsn",
    236: "tsl4",
    246: "clw",
    141: "snd",
    198: "src",
    235: "skt",
    39: "swvl1",
    40: "swvl2",
    41: "swvl2",
    138: "svo",
    155: "sd",
}

levelmap = {
    129: "surface",
    130: "model_level",
    131: "model_level",
    132: "model_level",
    133: "model_level",
}


codemap = {var: code for code, var in varmap.items()}
show_cdo = True


def _get_attrs(code):
    eccodes = ectable.table
    data = eccodes.loc[code]
    attrs = data.to_dict()
    attrs["code"] = code
    return attrs


def _get_row_by_date(code, date, df, frequency="hourly", era_id="E5", use_E1=True):
    """Find a row in an ERA5 catalog dataframe.

    The path to a file should be uniquely identified by code, date, frequency
    and era_id. If this selection is not unique, we assume, there are some files
    with the same content (which is the case with some data downloaded later).
    If the date is between 2000 and 2006, ERA5.1 date can be used (if available)
    identified by ``era_id=='E1'``.

    """
    date = pd.to_datetime(date)
    era_id_sel = era_id
    if use_E1 is True:
        if date >= pd.to_datetime("2000-01-01T00:00:00") and date <= pd.to_datetime(
            "2006-12-31T23:00:00"
        ):
            era_id_sel = "E1"
    # lt = "model_level"
    sel = (
        df.set_index(["code", "frequency", "time"]).sort_index().loc[(code, frequency)]
    )

    # try to select era5.1 if possible
    if not sel[sel.era_id == era_id_sel].empty:
        sel = sel[sel.era_id == era_id_sel]
    else:
        sel = sel[sel.era_id == era_id]

    if not sel.level_type.nunique() == 1:
        if code in levelmap:
            sel = sel[sel["level_type"] == levelmap[code]]
        else:
            raise Exception(
                "non unique selection: {}, {}".format(code, sel.level_type.unique())
            )

    # time index should be unique now...
    if not sel.index.is_unique:
        sel = sel[~sel.index.duplicated(keep="first")]
    sel = sel.sort_index()
    # select date index
    ix = sel.index.get_indexer([date], method="pad")
    return sel.iloc[ix]


def _get_file_by_date(code, date, df):
    """returns a filepath depending on code and data

    This subroutine is used to find the ERA5 grib file
    that contains the code and date.

    """
    return _get_row_by_date(code, date, df).path.values[0]


def _rename_variable(ds, codes):
    """change name and attributes of dataarray

    Trys to identifiy, rename and update attributes
    of an ERA5 variable in an xarray dataset. Since we don't
    exactly know how cdo converts the meta data and attributes
    of an ecmwf grib file, we try several ways to identify
    a variable, e.g, by name, code, etc...

    """
    if not isinstance(codes, list):
        codes = [codes]
    rename = {}
    for code in codes:
        attrs = _get_attrs(code)
        for var, da in ds.items():
            if attrs["short_name"] == var:
                rename[var] = varmap[code]
            elif "code" in da.attrs:
                rename[var] = varmap[da.code]
            if var in rename:
                da.attrs = _get_attrs(code)
    return ds.rename(rename)


def _open_catalog(url="/pool/data/Catalogs/mistral-era5.json"):
    import intake

    cat = intake.open_esm_datastore(url)
    level_types = ["model_level", "surface"]
    return cat.search(
        code=list(codemap.values()),
        level_type=level_types,
        dataType="an",
        frequency="hourly",
    )


def _get_catalog_df(url="/pool/data/Catalogs/mistral-era5.json"):
    df = _open_catalog(url).df
    df["time"] = pd.to_datetime(df.validation_date)
    df = df.set_index(["time"])
    return df


def _get_code(ident):
    if type(ident) == int:
        return ident
    elif ident in codemap:
        return codemap[ident]
    else:
        message = (
            "Unkown code or cf variable name: "
            + str(ident)
            + "   I know: "
            + str(codemap)
        )
        raise Exception(message)


def _clean_coords(ds, vcs=None):
    if vcs is None:
        vcs = ["hyai", "hybi", "hyam", "hybm"]
    for vc in vcs:
        ds[vc] = ds[vc].isel(time=0).drop("time").squeeze()
    ds["akgm"] = ds.hyai
    ds["bkgm"] = ds.hybi
    return ds


class ERA5:
    """
    Class for cmorizing original ERA5 GRIB data.

    Notes
    -----
    The cmorizer class only works with the intake catalog provided by DKRZ.

    """

    def __init__(
        self,
        catalog_url="/pool/data/Catalogs/mistral-era5.json",
        df=None,
        scratch=None,
        show_cdo=True,
    ):
        """

        Parameters
        ----------
        catalog_url : str
            Url of the ERA5 intake catalog.
        scratch : str
            Scratch directory for temporary cdo files.

        """
        self.catalog_url = catalog_url
        self.show_cdo = show_cdo
        self.scratch = init_tempdir(scratch)
        self.df = df
        if self.df is None:
            self.df = _open_catalog(catalog_url).df
        self.df["time"] = pd.to_datetime(self.df.validation_date)

    def _to_dataarray(
        self,
        codes,
        dates,
        df,
        parallel=False,
        cf_meta=True,
        use_cftime=True,
        chunks={"time": 1},
        **kwargs
    ):
        timesteps = self._get_timesteps(codes, dates, delayed=parallel)
        if parallel:
            import dask

            timesteps_ = dask.compute(timesteps)[0]
        else:
            timesteps_ = timesteps
        print(timesteps_)
        dsets = {
            code: xr.open_mfdataset(
                files, use_cftime=use_cftime, chunks=chunks, **kwargs
            )
            for code, files in timesteps_.items()
        }
        if cf_meta is True:
            for code, ds in dsets.items():
                dsets[code] = _rename_variable(ds, code)
        return xr.merge(dsets.values(), compat="override", join="override")

    def to_xarray(self, idents, dates, parallel=False, cf_meta=True):
        """Create an xarray dataset from ERA5 GRIB data at DKRZ.

        Parameters
        ----------
        idents : str or int
            Code or variable name (cf name) that should be converted.
        dates : date or list of dates in ISO 8601 format.
            A single date or list of dates for which the variables should
            be converted.
        parallel : Bool
            Use dask delayed for converting several variables and dates.
        cf_meta : bool
            Rename variables to CF standard names.

        Returns
        -------
        Dataarray

        """
        if not isinstance(idents, list):
            idents = [idents]
        idents = [_get_code(ident) for ident in idents]
        ds = self._to_dataarray(
            idents, dates, self.df, parallel=parallel, cf_meta=cf_meta
        )
        return ds

    # own cdo call here due to https://github.com/Try2Code/cdo-bindings/issues/34
    def _cdo_call(self, options="", op="", input="", output="temp", print_command=True):
        if output is None:
            output = ""
        elif output == "temp":
            output = tempfile.TemporaryDirectory(dir=self.scratch).name
            tempfiles.append(output)
        if isinstance(input, list):
            input = " ".join(input)
        call = "{} {} {} {} {}".format(cdo_exe, options, op, input, output)
        if print_command is True:
            print(call)
        stdout = subprocess.Popen(
            call, shell=True, stdout=subprocess.PIPE
        ).stdout.read()
        if output:
            return output
        return stdout

    def _seldate(self, code, date):
        """select a single timestep of ERA5 data

        Returns the path to a temporay file that contains only
        a single timestep of ERA5 data.

        """
        filename = _get_file_by_date(code, date, self.df)
        return self._cdo_call(op="seldate,{}".format(date), input=filename)

    def _to_regular(self, filename, table=None):
        """converts ecmwf spectral grib data to regular gaussian netcdf.

        cdo is used to convert ecmwf grid data to netcdf depending on the gridtype:
        For 'gaussian_reduced': cdo -R
            'spectral'        : cdo sp2gpl

        This follows the recommendation from the ECMWF Era5 Documentation.
        We also invert the latitudes to stick with cmor standard.

        """
        if table is None:
            table = ""
        # from cdo import Cdo
        # cdo = Cdo(tempdir=scratch)
        gridtype = self._gridtype(filename)
        if gridtype == "gaussian_reduced":
            # return self.cdo.copy(options='-R -f nc', input=filename)
            # gaussian = cdo.setgridtype('regular', options='-f nc', input=filename)
            # gaussian = _cdo_call(op='set')
            gaussian = self._cdo_call(
                op="setgridtype,regular", options="-f nc " + table, input=filename
            )
        elif gridtype == "spectral":
            # gaussian = Cdo(tempdir=scratch).sp2gpl(options='-f nc', input=filename)
            # gaussian = _convert_with_cdo(filename, op='sp2gpl')
            gaussian = self._cdo_call(
                op="sp2gpl", options="-f nc " + table, input=filename
            )
        elif gridtype == "gaussian":
            # gaussian =  cdo.copy(options='-f nc', input=filename)
            gaussian = self._cdo_call(
                op="copy", options="-f nc " + table, input=filename
            )
        else:
            raise Exception(
                "unknown grid type for conversion to regular grid: {}".format(gridtype)
            )
        return self._cdo_call(op="invertlat", input=gaussian)
        # return gaussian

    def _split_time(self, filename, scratch=None):
        from cdo import Cdo

        base = os.path.basename(filename)
        if scratch is None:
            scratch = os.path.join(os.environ["SCRATCH"], ".cdo-scratch")
            if not os.path.exists(scratch):
                os.mkdir(scratch)
        # output=os.path.join(scratch, base)
        return Cdo().splitsel(
            "1,0,5", input=filename, output=os.path.join(scratch, base)
        )
        # _cdo_call(op='splitsel,1', input=filename, output=output)
        # return output

    def _get_timestep(self, code, date, gaussian=True):
        f = self._seldate(code, date)
        if gaussian is True:
            return self._to_regular(f)
        return f

    def _extract_timesteps(self, codes, dates, gaussian=True, delayed=False):
        from itertools import product

        if not isinstance(dates, list):
            dates = [dates]
        if not isinstance(codes, list):
            codes = [codes]
        timesteps = {code: [] for code in codes}
        if delayed is True:
            from dask import delayed
        else:

            def delayed(x):
                return x

        for code, date in product(codes, dates):
            reg = delayed(self._get_timestep)(code, date, gaussian)
            timesteps[code].append(reg)
        return timesteps

    def _get_timesteps(self, idents, dates, gaussian=True, delayed=False):
        """Returns a list of cmorized files"""
        idents = [_get_code(ident) for ident in idents]
        return self._extract_timesteps(
            idents, dates, gaussian=gaussian, delayed=delayed
        )

    def wind(self, dates, parallel=False, cf_meta=True):
        """Create an xarray dataset of wind variables from ERA5 GRIB data at DKRZ.

        Parameters
        ----------
        dates : date or list of dates in ISO 8601 format.
            A single date or list of dates for which the variables should
            be converted.
        parallel : Bool
            Use dask delayed for converting several variables and dates.
        cf_meta : bool
            Rename variables to CF standard names.

        Returns
        -------
        Dataarray

        """
        return self._wind(dates, parallel=parallel, cf_meta=cf_meta)

    def _atmosphere(self, dates, parallel=False, cf_meta=True):
        variables = ["ta", "hus", "ps", "tos", "sic", "clw", "snd"]
        if parallel is True:
            from dask import delayed
        else:

            def delayed(x):
                return x

        ds = delayed(self.to_xarray)(
            variables, dates, parallel=parallel, cf_meta=cf_meta
        )
        wind = delayed(self.wind)(dates, parallel=parallel, cf_meta=cf_meta)
        if parallel is True:
            import dask

            res = dask.compute(ds, wind)
        else:
            res = (ds, wind)
        # print(ds)
        # print(wind)
        return xr.merge(res, join="override", compat="override")

    def _dynamics(self, dates, parallel=False, cf_meta=True):
        if not isinstance(dates, list):
            dates = [dates]
        if parallel is True:
            from dask import delayed
        else:

            def delayed(x):
                return x

        res = []
        for date in tqdm(dates):
            res.append(delayed(self._atmosphere)(date))
        if parallel is True:
            import dask

            res_ = dask.compute(res)[0]
        else:
            res_ = res
        return xr.concat(res_, dim="time")

    def _fx(self, date=None, cf_meta=True):
        """static variables

        Static variables are supposed to be time dependent.

        """
        variables = ["orog", "sftlf"]
        if date is None:
            date = "1979-01-01T00:00:00"
        return self.to_xarray(variables, date, cf_meta=cf_meta).drop("time").squeeze()

    def gfile(
        self,
        dates,
        parallel=False,
        cf_meta=True,
        clean_coords=True,
        add_fx=True,
        fx_date=None,
    ):
        """Create an ERA5 gfile dataset.

        Main function to convert ERA5 grib data to a regular gaussian Dataset
        containing all variables required for REMO preprocessing.

        Parameters
        ----------
        dates : date or list of dates in ISO 8601 format.
            A single date or list of dates for which the variables should
            be converted.
        parallel : Bool
            Use dask delayed for converting several variables and dates.
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
        gds = self._dynamics(dates, parallel=parallel, cf_meta=cf_meta)
        if add_fx is True:
            gds = xr.merge(
                [gds, self._fx(date=fx_date, cf_meta=cf_meta)],
                join="override",
                compat="override",
            )
        if clean_coords is True:
            gds = _clean_coords(gds)
        return gds

    def _convert_with_cdo(self, f, op):
        """Convert a single file into NetCDF format."""
        file = os.path.basename(f)
        path = os.path.dirname(f)
        output = os.path.join(path, file + ".nc")
        self._cdo_call(options="-f nc", op=op, input=f, output=output)
        return output

    def _griddes(self, filename):
        griddes = (
            self._cdo_call(options="", op="griddes", input=filename, output=None)
            .decode("utf-8")
            .split("\n")
        )
        griddes = {
            entry.split("=")[0].strip(): entry.split("=")[1].strip()
            for entry in griddes
            if "=" in entry
        }
        return griddes  # ['gridtype']

    def _gridtype(self, filename):
        return self._griddes(filename)["gridtype"]

    def _wind_by_date(self, date):
        """compute wind from vorticity and divergence"""
        vorticity = "svo"
        divergence = "sd"
        vort_tmp = self._get_timestep(codemap[vorticity], date, gaussian=False)
        div_tmp = self._get_timestep(codemap[divergence], date, gaussian=False)
        # ds = xr.merge([vort_ds, div_ds])
        merge = self._cdo_call(op="merge", input=[vort_tmp, div_tmp])
        uv = self._cdo_call(op="dv2uvl", options="-f nc", input=merge)
        uv = self._cdo_call(op="invertlat", input=uv)
        return uv

    def _wind_timesteps(self, dates, delayed=False):
        if not isinstance(dates, list):
            dates = [dates]
        uv_dates = []
        if delayed is True:
            from dask import delayed
        else:

            def delayed(x):
                return x

        for date in dates:
            uv_date = delayed(self._wind_by_date)(date)
            uv_dates.append(uv_date)
        return uv_dates

    def _wind(self, dates, parallel=False, cf_meta=True):
        files = self._wind_timesteps(dates, delayed=parallel)
        if parallel:
            import dask

            files_ = dask.compute(files)[0]
        else:
            files_ = files
        ds = xr.open_mfdataset(
            files_,
            chunks={"time": 1},
            coords="minimal",
            use_cftime=True,
            compat="override",
            data_vars="minimal",
        )
        if cf_meta is True:
            return _rename_variable(ds, [codemap[var] for var in ["ua", "va"]])
        return ds
