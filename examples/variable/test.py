

import xarray as xr

import pint_xarray

ds = xr.Dataset(
          {
            "a": (("lon", "lat"), [[11.84, 3.12, 9.7], [7.8, 9.3, 14.72]]),
               "b": (("lon", "lat"), [[13, 2, 7], [5, 4, 9]], {"units": "K"}),
              },
             coords={"lat": [10, 20, 30], "lon": [74, 76]},
           )

print(ds.b)
b = ds.b.pint.quantify()

#print(b)
