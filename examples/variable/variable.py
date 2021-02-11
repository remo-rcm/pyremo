

from pyremo import remo_ds

filename = '../../tests/data/e056111t2006010100'

ds = remo_ds.open_remo_dataset(filename, update_meta=True)

print(ds)
