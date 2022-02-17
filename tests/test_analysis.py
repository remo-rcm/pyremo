import pytest
import pyremo as pr
import xarray as xr
import numpy as np
import cartopy.crs as ccrs

def test_analysis():

    data_tas = '/work/ch0636/g300100/remo_test_data/tas/tas_e056000m200601-200612.nc'
    data_pr = '/work/ch0636/g300100/remo_test_data/pr/pr_e056000m200601-200612.nc'
    data_lsm = '/work/ch0636/g300100/remo_test_data/lsm/lsm_e056000t1979010106.nc'

    obs_tas = '/work/ch0636/g300100/eobs_test_data/tas/tg_ens_mean_0.25deg_reg_v22.0e_200601-200612.nc'
    obs_pr= '/work/ch0636/g300100/eobs_test_data/pr/rr_ens_mean_0.25deg_reg_v22.0e_200601-200612.nc'

    remo_bla = pr.analysis.open_mfdataset(data_lsm)
    remo_temp2 = pr.analysis.remo_preprocess(pr.analysis.open_mfdataset(data_tas))
    remo_apr = pr.analysis.remo_preprocess(pr.analysis.open_mfdataset(data_pr))

    remo_lsm = xr.where(remo_bla.BLA.isel(time=0) > 0, 1, 0)
    remo_temp2 = (remo_temp2.TEMP2-273.15)
    #remo_tas.attrs['units'] = 'Celsius'
    remo_apr = remo_apr.APRL*24
    remo_tas = xr.merge([remo_temp2, remo_lsm]).rename({'TEMP2':'tas', 'BLA':'mask'})
    remo_pr = xr.merge ([remo_apr, remo_lsm]).rename({'APRL':'pr', 'BLA':'mask'})

    eobs_tg = pr.analysis.open_mfdataset(obs_tas)
    eobs_rr = pr.analysis.open_mfdataset(obs_pr)
    eobs_lsm = xr.where(~np.isnan(eobs_tg.isel(time=0)),1,0).rename({'tg':'mask'})
    eobs_tas = xr.merge([eobs_tg, eobs_lsm]).rename({'tg':'tas'})
    eobs_pr = xr.merge([eobs_rr, eobs_lsm]).rename({'rr':'pr'})

    pole = (remo_bla.rotated_latitude_longitude.grid_north_pole_longitude,  
            remo_bla.rotated_latitude_longitude.grid_north_pole_latitude)
    crs =  ccrs.RotatedPole(*pole)

    plotting_settings = {"transform" : ccrs.PlateCarree(),
                         "extent" : {"extents" : [remo_bla.rlon.min(), remo_bla.rlon.max(), 
                                                  remo_bla.rlat.min(), remo_bla.rlat.max()],
                                     "crs" : crs},
                         "projection" : crs,
                         "borders" : True,
                         "xlocs" : range(-180,180,10),
                         "ylocs" : range(-90,90,10),
                         "figsize" : (14,10),
                         "aspect" : "auto"}

    compare_eobs_tas = pr.analysis.compare_seasons(remo_tas, eobs_tas).compute()
    
    eobs_tas_season_plot = pr.analysis.plot_seasons(compare_eobs_tas.tas,
                                                    **plotting_settings)

    compare_eobs_pr = pr.analysis.compare_seasons(remo_pr, eobs_pr).compute()
    eobs_pr_season_plot = pr.analysis.plot_seasons(compare_eobs_pr.pr, 
                                                   **plotting_settings,
                                                   vmin=-2, vmax=2, cmap='bwr_r')

