"""
This file is responsible for converting UniProt IDs to Entrez IDs using BioDBNet

TODO: Convert this into a class for easier use
    - Maybe user will input a directory containing all sqt files?
TODO: Enable multiprocessing
    - Use process_sqt_file() function as a wrapper

TODO: Determine how user will input data
    - This may be a "plugin" for main_proteomics.py

"""
from bioservices import BioDBNet
import csv
import os
from pathlib import Path
import sys
import tqdm

from pprint import pprint

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


def parse_uniprot(uniprot_id: str) -> str:
    """
    This function is responsible for collecting the first-index field from a pipe-separated string

    EXAMPLE:
        sp|Q16832|DDR2_HUMAN -> Q16832
        sp|Q92736|RYR2_HUMAN -> Q92736
        tr|H0YL77|H0YL77_HUMAN -> H0YL77
    """
    # Split the string on the pipe
    uniprot_id = uniprot_id.split("|")[1]

    # Return the first index
    return uniprot_id


def read_uniprot_ids(input_csv_file: Path) -> dict:
    """
    This function is responsible for collecting the UniProt IDs from the input csv file

    Documentation: https://crux.ms/file-formats/sqt-format.html

    We must perform the following
    1. Skip every line starting with an "H"
    2. Skip every line starting with an "M"
    3. Collect every line starting with an "S". This contains the ion intensity
        - Collect the next line that starts with an "L". This contains the UniProt ID
    4. Repeat steps 2 and 3 until the end of the file
    """
    # Create an empty dictionary to store the UniProt IDs and ion intensities
    # uniprot_ids: dict[str, list[list[str]]
    # ion_intensities: dict[str, list[str]]
    uniprot_data: dict = {
        "ion_intensities": [],
        "uniprot_ids": [],
    }

    # Open the file
    with open(input_csv_file, "r") as i_stream:

        # Iterate through the lines in the file
        for line in i_stream:
            # If the line starts with an "S", collect it
            if line.startswith("S"):
                # Collect the 8th tab-separated field (index of 7)
                # Save this in the uniprot_ids dictionary
                uniprot_data["ion_intensities"].append(line.split("\t")[7])

                # Find the next "locus" line that does not start with "decoy_"
                # Save this in the uniprot_ids dictionary
                locus_found = False
                while not locus_found:
                    try:
                        line = next(i_stream)
                    except StopIteration:
                        break
                    # Only accept lines starting with "L" (locus)
                    if line[0] == "L":
                        # If the line starts with "decoy_", skip it
                        uniprot_id = line.split("\t")[1]
                        if uniprot_id.startswith("decoy_"):
                            continue  # Goes to the next line
                        else:
                            # Get the UniProt ID from the line
                            uniprot_data["uniprot_ids"].append(
                                parse_uniprot(uniprot_id)
                            )
                            locus_found = True  # Exits the while loop

    return uniprot_data


def convert_ids(uniprot_data: dict) -> dict:
    """
    This function is responsible for converting a list of uniprot IDs to Gene IDs
    """
    # Create a new BioDBNet object
    biodbnet = BioDBNet()
    # Create an empty list to store the gene IDs
    uniprot_data["gene_ids"] = []

    for i in tqdm.tqdm(range(0, len(uniprot_data["uniprot_ids"]), 500)):
        # for i in range(0, len(uniprot_data["uniprot_ids"]), 500):
        # Get the next 500 UniProt IDs from uniprot_data
        id_lookup = uniprot_data["uniprot_ids"][i : i + 500]

        # Convert the UniProt IDs to gene IDs
        gene_ids = biodbnet.db2db(
            "UniProt Accession", "Gene ID", input_values=id_lookup
        )

        # Add the gene IDs to the list
        uniprot_data["gene_ids"].extend(gene_ids["Gene ID"])

    # Return the updated dictionary
    return uniprot_data


def write_data(uniprot_data: dict, output_file_path: Path) -> None:
    """
    This function is responsible for writing a dictionary to a csv file

    The dictionary has the following keys:
    1. uniprot_ids
    2. gene_ids
    3. ion_intensities

    Each key has a list of values
    """
    # Open the output file
    with open(output_file_path, "w") as o_stream:
        # Create a csv writer
        writer = csv.writer(o_stream, delimiter=",")

        # Write the header
        writer.writerow(["uniprot_id", "gene_id", "ion_intensity"])

        # Write the data
        for uniprot_id, gene_id, ion_intensity in zip(
            uniprot_data["uniprot_ids"],
            uniprot_data["gene_id"],
            uniprot_data["ion_intensities"],
        ):
            writer.writerow([uniprot_id, gene_id, ion_intensity])


def find_sqt_file(input_dir: Path) -> set[str]:
    """
    This function is responsible for recursively finding the sqt file in the input directory
    It will return a set of all the sqt files in the input directory
    """
    # Create a list to store the files
    sqt_files: set[str] = set()

    # Recursively loop through the files in the input directory
    for root, dirs, files in os.walk(input_dir):

        # Loop through the files
        for file in files:

            # If the file is a sqt file, add it to the list
            if file.endswith(".sqt"):
                sqt_files.add(os.path.join(root, file))

    return sqt_files


def process_sqt_file(input_sqt_file: Path) -> dict:
    """
    This function is responsible for processing the sqt file
    """
    uniprot_ids: dict = read_uniprot_ids(input_sqt_file)

    # Convert uniprot IDs to a new list of gene IDs
    gene_ids: dict = convert_ids(uniprot_ids)

    # Return gene ids so they can be written to a file
    return gene_ids


def main(args: list[str]):
    """
    This is the main function for the script
    """
    config = project.configs
    input_dir: Path = Path(config.configdir, "uniprot_gene_id_conversion")
    output_file: Path = Path(config.outputdir, "gene_id_intensity_map.csv")

    # Get a list of sqt files
    sqt_files: set[str] = find_sqt_file(input_dir)
    gene_ids: dict = {}
    for i, file in enumerate(sqt_files):
        file_path: Path = Path(file)
        print(f"Processing {file} ({i + 1} / {len(sqt_files)})")

        # Extend the gene_ids dictionary with the new data
        gene_ids.update(process_sqt_file(file_path))

        print("")

    # Write the data to a csv file
    write_data(gene_ids, output_file)


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
