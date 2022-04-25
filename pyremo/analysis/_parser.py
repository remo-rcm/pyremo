import argparse


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Select path to REMO output data.")
    parser.add_argument(
        "-op",
        "--output_path",
        dest="output_path",
        nargs="?",
        default=".",
        help="Select output plotting destination. Default is the current directory.",
    )
    # parser.add_argument('-obs',
    #                    '--observational_data',
    #                    nargs='+',
    #                    default=['HYRAS', 'CRU_TS4', 'EOBS'],
    #                    help='Select observational data fro comparision. Choose between HYRAS, CRU_TS4 and EOBS in any combination. Default are all three datasets.')
    parser.add_argument(
        "-tr",
        "--time_range",
        dest="time_range",
        nargs="+",
        default=["1980", "2010"],
        help="Select analysis time range. Default is from 1980 to 2010.",
    )
    return parser


parser = create_parser()
args = parser.parse_args()
