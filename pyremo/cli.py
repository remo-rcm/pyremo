"""Console script for pyremo."""
import argparse
import sys
from .prsint import _defaults as dflt


def create_parser():
    """Creates parser for command line tool."""
    return argparse.ArgumentParser()


def main():
    """Console script for pyremo."""
    parser = argparse.ArgumentParser()
    parser.add_argument("_", nargs="*")
    args = parser.parse_args()

    print("Arguments: " + str(args._))
    print("Replace this message by putting your code into pyremo.cli.main")
    return 0


def prsint_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", metavar="input file", nargs="+")
    parser.add_argument(
        "-v",
        "--variables",
        dest="variables",
        nargs="+",
        help="list of variables to interpolate (default = {})".format(dflt.variables),
        default=dflt.variables,
    )
    parser.add_argument(
        "-p",
        "--plevs",
        dest="plevs",
        nargs="+",
        type=int,
        help="list of pressure levels to interpolate to (default = {})".format(
            dflt.plevs
        ),
        default=dflt.plevs,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="output format",
        choices=["plevs", "input"],
        default="plevs",
    )
    parser.add_argument(
        "-id",
        "--id",
        dest="id",
        help="experiment id for output file naming",
        default="000000",
    )
    parser.add_argument(
        "-cdo",
        "--cdo_options",
        dest="cdo_options",
        help="options for using cdo to read input",
        default="",
    )
    return parser


def prsint():
    """Console script for prsint."""
    parser = prsint_parser()
    args = parser.parse_args()

    print("Arguments: " + str(args))
    # druint.test()
    # prsint(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
