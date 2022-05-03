# flake8: noqa
#!/usr/bin/env python
# SBATCH --job-name=double-nesting
# SBATCH --partition=shared
# SBATCH --nodes=1
# SBATCH --ntasks=24
# SBATCH --time=24:00:00
# SBATCH --mem-per-cpu=2600
# SBATCH --account=ch0636
# SBATCH --nice=1000


import glob

import dask
import xarray as xr
from dask.distributed import Client, progress

import pyremo as pr
from pyremo.preproc import double_nesting as dn


def process_file(file, em, hm, vc, surflib, path=None, write=True):
    ds = pr.preprocess(xr.open_dataset(file))
    dn.init(ds, em, hm, vc, surflib)
    dsa = dn.remap_timestep(ds)
    dn.deallocate()
    if write is True:
        return dn.write_timestep(dsa, path)
    return dsa


def process_files(files, em, hm, vc, surflib, path=None, write=True):
    results = []
    for f in files:
        results.append(dask.delayed(process_file)(f, em, hm, vc, surflib, path, write))
    return results


def main():
    surflib = pr.preprocess(
        xr.open_dataset(
            "/work/ch0636/g300046/preproc_test/doublenest/lib/lib_EUC-0275_frac.nc"
        )
    ).squeeze(drop=True)

    path = "/scratch/g/g300046/197901/xa"
    print("writing to: {}".format(path))

    em = pr.domain_info("EUR-11")
    hm = pr.domain_info("EUC-0275")
    # hm_ds = pr.remo_domain("EUC-0275")  # , dummy='topo')
    # em_ds = pr.remo_domain("EUR-11")  # , dummy='topo')
    vc = pr.vc.tables["vc_27lev"]

    data_dir = "/mnt/lustre01/scratch/g/g300046/197901/e056000t197901*"

    files = glob.glob(data_dir)
    files.sort()

    return process_files(files, em, hm, vc, surflib, path)


if __name__ == "__main__":
    # client = Client(threads_per_worker=1, n_workers=24, dashboard_address=None)
    print("starting client...")
    with Client(threads_per_worker=1, n_workers=24, dashboard_address=None) as client:
        print(client)
        results = main()
        results_ = dask.persist(*results)
        progress(*results_)
        results_ = dask.compute(*results_)
        # client.shutdown()
