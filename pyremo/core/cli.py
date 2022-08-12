import argparse

soil_default = [
    "TSL",
    "TSN",
    "TD3",
    "TD4",
    "TD5",
    "GLAC",
    "TD",
    "TDCL",
    "WS",
    "WL",
    "SN",
]


def replace_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "target", metavar="target file to which the soil should be added"
    )
    parser.add_argument(
        "source", metavar="source file from which the soil should be copied"
    )
    parser.add_argument(
        "-v",
        "--variables",
        dest="variables",
        nargs="+",
        help="list of variables to replace (default = {})".format(soil_default),
        default=soil_default,
    )
    return parser


def replace_variables(args):
    pass
