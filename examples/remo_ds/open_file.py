

from pyremo import remo_ds as rds


filename = "../../tests/data/e056111t2006010100"

ds = rds.open_remo_dataset(filename, update_meta=True, parse_dates=True)
