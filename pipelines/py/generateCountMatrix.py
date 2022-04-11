#!/usr/bin/python3
import argparse
import sys
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage

# from rpy2.robjects import r, pandas2ri
# import rpy2.robjects as ro
# from rpy2.robjects.conversion import localconverter
# from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# pandas2ri.activate()

# limma = importr("limma")
tidyverse = importr("tidyverse", lib_loc="/home/jupyteruser/rlibs")
# edgeR = importr("edgeR")
# genefilter = importr("genefilter")
# biomaRt = importr("biomaRt")
# sjmisc = importr("sjmisc")

# automatically convert ryp2 dataframe to Pandas dataframe
string = """
c
"""
genCountMatrixio = SignatureTranslatedAnonymousPackage(string, "genCountMatrixio")


def main(argv):
    # TODO: Fix description
    parser = argparse.ArgumentParser(
        prog="generateCountMatrix.py",
        description="This description must be fixed",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-d",
        "--directory_input",
        type=str,
        required=True,
        dest="directory_input",
        help="The input directory to get files from",
    )
    parser.add_argument(
        "-o",
        "--directory_output",
        type=str,
        required=True,
        dest="directory_output",
        help="The directory to save outputs to",
    )
    parser.add_argument(
        "-t",
        "--technique",
        type=str,
        required=True,
        dest="technique",
        help="The experiment technique",
    )
    args = parser.parse_args()
    input_dir = args.directory_input
    output_dir = args.directory_output
    technique = args.technique

    # TODO: Remove this after verifying argparse works properly
    # try:
    #     opts, args = getopt.getopt(
    #         argv, "hi:o:t:", ["input_dir=", "output_dir=", "technique="]
    #     )
    # except getopt.GetoptError:
    #     print(
    #         "python3 proteomics_gen.py -i <input directory> -o <output_directory> -t <technique>"
    #     )
    #     sys.exit(2)
    # for opt, arg in opts:
    #     if opt == "-h":
    #         print(
    #             "python3 proteomics_gen.py -i <input_directory> -o <output_directory> -t <technique>"
    #         )
    #         sys.exit()
    #     elif opt in ("-i", "--input_dir"):
    #         input_dir = arg
    #     elif opt in ("-o", "--output_dir"):
    #         output_dir = arg
    #     elif opt in ("-t", "--technique"):
    #         technique = arg

    print(f"Input directory is '{input_dir}'")
    print(f"Output directory is '{output_dir}'")
    print(f"Active gene determination technique is '{technique}'")

    genCountMatrixio.genCountMatrix_main(input_dir, output_dir, technique)
    print("Out")


if __name__ == "__main__":
    main(sys.argv)
