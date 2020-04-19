#! /usr/bin/env python
# -*- coding: utf-8 -*-


#################################
#
#   MAIN
#
#################################
if __name__ == "__main__":

    import argparse
    import numpy as np
    from ztflc import forcephotometry

    # ================= #
    #   Options         #
    # ================= #
    parser = argparse.ArgumentParser(
        description=""" run the interactive plotting of a given cube""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("infile", type=str, default="None", help="Target Name")

    parser.add_argument(
        "--download",
        action="store_true",
        default=False,
        help="Download the data if necessary ?",
    )

    parser.add_argument(
        "--ndlprocess",
        type=int,
        default=4,
        help="Number of parallel processing for the download",
    )

    # // Other stuffs
    parser.add_argument(
        "--quiet", action="store_true", default=False, help="Set the verbose to False"
    )

    # Parsing
    args = parser.parse_args()

    # ================= #
    #   Script          #
    # ================= #
    forcephotometry.run_forcephotometry(
        args.infile,
        download_files=args.download,
        nprocess=args.ndlprocess,
        verbose=~args.quiet,
    )
