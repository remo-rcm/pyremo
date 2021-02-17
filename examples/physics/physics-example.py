
from pyremo import remo_ds
from pyremo.physics import compute_arfgm, relative_humidity, specific_humidity, liquid_water_content, pressure_at_model_levels

import pyintorg

import intorg

filename = '../../tests/data/e056111t2006010100'
ds = remo_ds.open_remo_dataset(filename, update_meta=True)

#print(ds)

ak = ds.hyai
bk = ds.hybi

akh = 0.5 * (ak[:-1] + ak[1:])
bkh = 0.5 * (bk[:-1] + bk[1:])

print(10*'-')
print(ak)
print(bk)
print(10*'-')
print(akh)
print(bkh)
print(10*'-')
print(ds.hyam)
print(ds.hybm)
print(10*'-')

#relhum = relative_humidity(ds.T[0].T, ds.QD[0].T, ds.PS[0].T, ak, bk)
relhum = compute_arfgm(ds.T[0].T, ds.QD[0].T, ds.PS[0].T, ds.hyam, ds.hybm)
relhum_f = pyintorg.relative_humidity(ds.QD[0].T, ds.T[0].T, ds.PS[0].T, ds.hyai.values, ds.hybi.values, qwgm=None)

ph = pressure_at_model_levels(ds.PS[0].T, ds.hyam, ds.hybm)
ph_f = intorg.mo_test.compute_pressure(ds.PS[0].T, ak.values, bk.values)
#ph_f = intorg.mo_test.compute_pressure(ds.PS[0].T, ak, bk)

print(ph.shape)
print(ph[:,:,26])
print(ph_f.shape)
print(ph_f[:,:,26])

print(ds.PS[0].T.shape)
print(ds.PS[0].T.values)

print(ph == ph_f)

print(relhum.shape)
print(relhum[:,:,0])
print(relhum_f.shape)
print(relhum_f[:,:,0])
##
##print(relhum.max())
##print(relhum_f.max())
###
##
##import numpy as np
##print(np.isfortran(relhum))
##print(np.isfortran(relhum_f))
##
qw = liquid_water_content(ds.T[0].T, relhum, ds.PS[0].T, ds.hyam, ds.hybm)
qd = specific_humidity(ds.T[0].T, relhum, ds.PS[0].T, ds.hyam, ds.hybm)

print(10*'-')
print(qd[:,:,0])
print(ds.QD[0,0,:,:].T.values)
print(qw[:,:,0])
print(ds.QW[0,0,:,:].T.values)
print(10*'-')
