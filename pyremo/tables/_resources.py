from os import path as op

import pandas as pd
import pooch

base_url = "https://raw.githubusercontent.com/remo-rcm/tables/main/"

cache_url = "~/.pyremo"

DOMAIN_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "domains/",
    registry={
        "domains.csv": None,
    },
)


CODE_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "code-list/",
    registry={
        "code-list.csv": None,
        "flake.csv": None,
    },
)


VC_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "vc/",
    registry={
        "vc_27lev.csv": None,
        "vc_40lev.csv": None,
        "vc_49lev.csv": None,
        "vc_49lev_nh_pt2000.csv": None,
        "vc_60lev_ecmwf.csv": None,
        "vc_91lev_ecmwf.csv": None,
        "vc_101lev_100m_pt27713.csv": None,
        "vc_144lev_equal_pdiff_pt6400.csv": None,
    },
)


def fetch_remote_table(name, resource):
    """
    uses pooch to cache files
    """

    # the file will be downloaded automatically the first time this is run.
    return resource.fetch(name)


def read_remote_table(name, resource, **kwargs):
    fname = fetch_remote_table(name, resource)

    return pd.read_csv(fname, **kwargs)


def read_remo_domain_tables():
    resource = DOMAIN_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(table, resource, index_col="domain_id")
        for table in resource.registry.keys()
    }


# def read_remo_code_tables():
#     resource = CODE_RESOURCE
#     return {
#         table.split(".")[0]: read_remote_table(table, resource, index_col="code")
#         for table in resource.registry.keys()
#     }


def read_remo_code_tables():
    filename = pooch.retrieve(
        op.join(base_url, "code-list", "table.csv"), known_hash=None, path=cache_url
    )
    return {"codes": pd.read_csv(filename, index_col="code")}


def read_remo_vc_tables():
    resource = VC_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(
            table, resource, comment="#", index_col=None
        )
        for table in resource.registry.keys()
    }
