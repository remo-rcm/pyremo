import numpy as np
import xarray as xr

from ..archive import archive
from ..core import codes
from ..core.remo_ds import preprocess as remo_preprocess
from .obs import height_correction

soil_temps = ["TS", "TSL", "TEMP2", "TSN", "TD3", "TD4", "TD5"]


def open_mfdataset(
    files,
    use_cftime=True,
    parallel=True,
    data_vars="minimal",
    chunks={"time": 1},
    coords="minimal",
    compat="override",
    drop=None,
    **kwargs
):
    """optimized function for opening large cf datasets.

    based on https://github.com/pydata/xarray/issues/1385#issuecomment-561920115

    """

    def drop_all_coords(ds):
        return ds.reset_coords(drop=True)

    ds = xr.open_mfdataset(
        files,
        parallel=parallel,
        decode_times=False,
        combine="by_coords",
        preprocess=drop_all_coords,
        decode_cf=False,
        chunks=chunks,
        data_vars=data_vars,
        coords=coords,
        compat=compat,
        **kwargs
    )
    return xr.decode_cf(ds, use_cftime=use_cftime)


class VariableSet:
    def __init__(self, varnames, name=None, long_name=None):
        self.varnames = varnames

    def extract(self, ds):
        res = xr.merge([ds[var] for var in self.varnames if var in ds])
        res.attrs = ds.attrs
        return res

    def stack(self, ds):
        return stack_variables(ds, self.varnames)


soil = VariableSet(soil_temps, name="soil_temperature", long_name="soil temperature")


def stack_variables(ds, varnames, name=None, long_name=None):
    found = [ds[var] for var in varnames if var in ds]
    dim = xr.DataArray(data=[var.name for var in found], dims="var", name="var")
    stack = xr.concat(found, dim=dim)
    if name is not None:
        stack.name = name
    if long_name is not None:
        stack.attrs["long_name"] = long_name
    return stack


def weighted_field_mean(ds, lon="rlon", lat="rlat", weights=None):
    """
    function to compute area-weighted spatial means
    """
    if weights is None:
        weights = np.cos(np.deg2rad(ds[lat]))
    return ds.weighted(weights).mean(dim=(lon, lat))


def daily_sum(da):
    """
    Function to compute daily sums with a simple groupby approach
    """

    return da.groupby("time.day").sum(dim="time")


def monthly_sum(da):
    """
    Function to compute dailymonthly means with a simple groupby approach
    """

    return da.groupby("time.month").mean(dim="time")


def seasonal_mean(da):
    """
    Optimized function to calculate seasonal averages from time series of monthly means

    based on: https://xarray.pydata.org/en/stable/examples/monthly-means.html
    """

    # Get number od days for each month
    month_length = da.time.dt.days_in_month
    # Calculate the weights by grouping by 'time.season'.
    weights = (
        month_length.groupby("time.season") / month_length.groupby("time.season").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby("time.season").sum().values, np.ones(4))

    # Calculate the weighted average
    return (da * weights).groupby("time.season").sum(dim="time")


class RemoExperiment:
    def __init__(self, path, time_range=None, preprocess=True):
        self.archive = archive.RemoArchive(path, force_parse=True)
        self.first_tfile = self.archive.first_tfile()

        self.ds = self._init_dataset(preprocess, time_range)

    def _init_dataset(self, preprocess=True, time_range=None):
        monthly = open_mfdataset(self.archive.monthly(time_range=time_range))
        if preprocess is True:
            monthly = remo_preprocess(monthly)
        lsm = xr.open_dataset(self.first_tfile).BLA.squeeze(drop=True)
        orog = xr.open_dataset(self.first_tfile).FIB.squeeze(drop=True)
        ds = xr.merge([monthly, lsm, orog])
        ds["mask"] = xr.where(lsm > 0, 1, 0)
        return ds

    @property
    def mask(self):
        return self.ds.mask

    @property
    def variables(self):
        return list(self.monthly.data_vars)

    def dataset(self, variables=None):
        if variables is None:
            ds = self.ds
        else:
            if not isinstance(variables, list):
                variables = [variables]
            names = [codes.get_dict(id)["variable"] for id in variables]
            ds = xr.merge([self.ds[v] for v in names])
        ds.attrs = self.monthly.attrs
        return ds

    def soil(self):
        return stack_variables(
            self.ds, soil_temps, name="soil_temperatures", long_name="soil temperatures"
        )

    def spin_up(self):
        return weighted_field_mean(self.soil())


def get_regridder(finer, coarser, method="bilinear", **kwargs):
    """
    Function to regrid data bilinearly to a coarser grid
    """

    import xesmf as xe

    return xe.Regridder(finer, coarser, method=method, **kwargs)


def compare_seasons(
    ds1, ds2, regrid="ds1", do_height_correction=False, orog1=None, orog2=None
):
    """
    Function to compare seasonal means of two datasets.

    Paramters
    ---------
    ds1 : xarray.Dataset
        First variable data for comparision. Temporal resolution has to be less than monthly.
        ds1 is mainly model output data.
        ds2 is subtracted from ds1.
    ds2 : xarray.Dataset
        Second variable data for comparision. Temporal resolution has to be less than monthly.
        ds2 is mainly observational or reanalysis data.
        ds2 is subtracted from ds1.
    regrid : {"ds1", "ds2"}, optional
        Denotes the dataset to be bilinearly regridded. Specify the dataset with the finer spatial resolution:

        - "ds1": Regrid ds1 to ds2's grid with coarser spatial resolution.
        - "ds2": Regrid ds2 to ds1's grid with coarser spatial resolution.
    do_height_correction : bool, optional
        If ``do_height_correction=True``, do a height correction on ds1 using two orography files orog1 and orog2.
    orog1 : xarray.Dataset, optional
        Use only if ``do_height_correction=True``.
        Specify a orography file referring to ds1.
    orog2 : xarray.Dataset, optional
        Use only if ``do_height_correction=True``.
        Specify a orography file referring to ds2.

    Returns
    -------
    seasonal_comparision : xarray.Dataset
        Spatial mean differences of two datasets

    """
    ds1 = ds1.copy()
    ds2 = ds2.copy()
    ds1_seasmean = seasonal_mean(ds1)
    ds2_seasmean = seasonal_mean(ds2)
    if regrid == "ds1":
        regridder = get_regridder(ds1, ds2)
        print(regridder)
        ds1_seasmean = regridder(ds1_seasmean)
    elif regrid == "ds2":
        regridder = get_regridder(ds2, ds1)
        print(regridder)
        ds2_seasmean = regridder(ds2_seasmean)

    if do_height_correction is True:
        orog1 = regridder(orog1)
        ds1_seasmean += height_correction(orog1, orog2)
    return ds1_seasmean - ds2_seasmean
    # return xr.where(ds1_seasmean.mask, ds2_seasmean - ds1_seasmean, np.nan)
