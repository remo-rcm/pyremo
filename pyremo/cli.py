"""Console script for pyremo."""
import argparse
import sys


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


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
