"""remo tables urls
"""

_code = {
    "basic": "https://raw.githubusercontent.com/remo-rcm/tables/master/code-list/code-list.csv",
    "flake": (
        "https://raw.githubusercontent.com/remo-rcm/tables/master/code-list/flake.csv"
    ),
}


# domain_tables = {'CORDEX': 'cordex-domains.csv',
#        'AUX': 'aux-domains.csv',
#        'CORDEX-FPS': 'cordex-fps.csv'}

_domain = {
    "CORDEX": "https://raw.githubusercontent.com/remo-rcm/tables/master/domains/cordex-domains.csv",
    "AUX": "https://raw.githubusercontent.com/remo-rcm/tables/master/domains/aux-domains.csv",
    "CORDEX-FPS": "https://raw.githubusercontent.com/remo-rcm/tables/master/domains/cordex-fps.csv",
}


_vc = {
    "27_standard": (
        "https://raw.githubusercontent.com/remo-rcm/tables/master/vc/vc_27lev.csv"
    ),
    "40_standard": (
        "https://raw.githubusercontent.com/remo-rcm/tables/master/vc/vc_40lev.csv"
    ),
    "49_standard": (
        "https://raw.githubusercontent.com/remo-rcm/tables/master/vc/vc_49lev.csv"
    ),
}


_table_address_dict = {"domain": _domain, "code": _code, "vc": _vc}
