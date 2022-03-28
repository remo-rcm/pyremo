
import pyremo as pr
import cordex as cx

# variables that should have a mask with fill values
fillvars = ["TSW", "SEAICE", "TSI"]


def get_grid(domain_info):
    return cx.create_dataset(**domain_info)



def encoding(da, missval):
    if np.isnan(da.values).any():
        return {"_FillValue": missval}
    else:
        return {"_FillValue": None}


def update_attrs(ds):
    for var, da in ds.items():
        try:
            attrs = pr.codes.get_dict(var)
            da.attrs = {}
            da.attrs["name"] = attrs["variable"]
            da.attrs["code"] = attrs["code"]
            da.attrs["description"] = attrs["description"]
            da.attrs["units"] = attrs["units"]
            # da.attrs['layer'] = attrs['layer']
            da.attrs["grid_mapping"] = "rotated_latitude_longitude"
            da.attrs["coordinates"] = "lon lat"
        except:
            pass
    return ds



def to_netcdf(
    ads,
    path="",
    expid="000000",
    template=None,
    tempfiles=None,
    missval=1.0e20,
    **kwargs
):
    """write dataset to netcdf

    by default, each timestep goes into a separate output file

    """
    if not os.path.isdir(path):
        os.makedirs(path)
    expand_time = [var for var, da in ads.items() if "time" in da.dims]
    if template is None:
        template = "a{}a{}.nc"
    dates, datasets = zip(*ads.groupby("time"))
    paths = [os.path.join(path, get_filename(date, expid, template)) for date in dates]
    dsets = []
    # expand time dimension only for variables not coordinates.
    for ds in datasets:
        for var, da in ds.items():
            if var in expand_time:
                ds[var] = da.expand_dims("time")
            if var in fillvars:
                ds[var].encoding["_FillValue"] = missval
            else:
                ds[var].encoding["_FillValue"] = None
        dsets.append(ds)
    # dsets = [dset.expand_dims('time') for dset in datasets]
    writer = xr.save_mfdataset(dsets, paths, **kwargs)
    if tempfiles is not None:
        for f in tempfiles:
            os.remove(f)
    return paths


def to_tar(files, tar_file, mode="w"):
    import tarfile

    try:
        from tqdm.notebook import tqdm
    except:

        def tqdm(x):
            return x

    tf = tarfile.open(tar_file, mode=mode)
    for f in tqdm(files, desc="creating tarfile"):
        tf.add(f, arcname=os.path.basename(f), recursive=False)
    tf.close()
    return tar_file


def encode(ds, missval=1.0e20):
    for var, da in ds.items():
        if var in expand_time:
            ds[var] = da.expand_dims("time")
        if var in fillvars:
            ds[var].encoding["_FillValue"] = missval
        else:
            ds[var].encoding["_FillValue"] = None