import pandas as pd
import pooch

base_url = "https://raw.githubusercontent.com/remo-rcm/tables/master/"

cache_url = "~/.pyremo"

DOMAIN_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "domains/",
    registry={
        "cordex-domains.csv": "8599e88ea84d6f844b9c5a7250325d0d91053d129334c09a508ad30c9a6e6c42",
        "aux-domains.csv": "2089d3f8eb9ec3a1abb0eddea40522f0297f4e7c8c0cdf11a6f90167caa7e5b3",
        "cordex-fps.csv": "02390f5255ba9a7b03a30654cc9ffd8944a95d7a076a6e729c35600ebd1cb72d",
    },
)


CODE_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "code-list/",
    registry={
        "code-list.csv": "e252298d78714be9346c27cf2a9048c4ddb1a0397b5cbf273b934f6b495eab01",
        "flake.csv": "eac78921934110d52513c042994186614816aaae373a65383de420c3080ad876",
    },
)


VC_RESOURCE = pooch.create(
    # Use the default cache folder for the OS
    path=cache_url,  # pooch.os_cache("cordex"),
    # The remote data is on Github
    base_url=base_url + "vc/",
    registry={
        "vc_27lev.csv": "ec14492d3cd9449f852ce90edd8aa3e56f37672be5e7e10d7f8672cc84764cb9",
        "vc_40lev.csv": "b2f946ccb9b4fd9a349a16396660efbbb89904fa42b5a1aaab113193a7ee0bad",
        "vc_49lev.csv": "d81595c437328b6aeba86c7ed7b35b50b275e0b43f0284f93471f99a5e81e41b",
        "vc_49lev_nh_pt2000.csv": "c5b9969bb1836f06ee87d3f3f80c06a489ff19d74542684bb45490916932541b",
        "vc_60lev_ecmwf.csv": "95afc9873c60d973a980861420cd2effaa4df2a9f913c4d6a7bea1f944a6bd4e",
        "vc_91lev_ecmwf.csv": "8abe07a0a2f844e0aa2691d03110b8ba2bdf295e2397e93d8267060be8d7a274",
        "vc_101lev_100m_pt27713.csv": "b345c43cc86e58abe64aedd7d3e88543217c7c9f0ffb85fa166c6c105fb33301",
        "vc_144lev_equal_pdiff_pt6400.csv": "0e10ceb27cb047f88322f860e6dfd8f77a9b90e431adcf88277996a8d2cbe021",
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


def read_domain_table(name):

    return read_remote_table(name, resource=DOMAIN_RESOURCE, index_col="short_name")


def read_remo_domain_tables():

    resource = DOMAIN_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(table, resource, index_col="short_name")
        for table in resource.registry.keys()
    }


def read_remo_code_tables():

    resource = CODE_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(table, resource, index_col="code")
        for table in resource.registry.keys()
    }


def read_remo_vc_tables():

    resource = VC_RESOURCE
    return {
        table.split(".")[0]: read_remote_table(
            table, resource, comment="#", index_col=None
        )
        for table in resource.registry.keys()
    }
