def horizontal_dims(da):
    """Returns the names of the horizontal dimensions."""
    for dim in da.dims:
        if "lon" in dim:
            lon_dim = dim
        if "lat" in dim:
            lat_dim = dim
    return (lon_dim, lat_dim)


def encode(ds, coord_fill_value=None):
    """Encode a dataset for writing to NetCDF."""
    for coord in ds.coords:
        ds.coords[coord].encoding["_FillValue"] = coord_fill_value
    return ds


def read_yaml(filename):
    import yaml

    with open(filename, "r") as f:
        return yaml.safe_load(f)
