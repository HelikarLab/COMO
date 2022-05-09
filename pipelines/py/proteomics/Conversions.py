"""
This file is responsible for converting UniProt IDs to Entrez IDs using BioDBNet
"""
from bioservices import BioDBNet
import csv
import os

from pprint import pprint


def read_uniprot_ids(input_csv_file: str) -> dict:
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
        "uniprot_ids": [],
        "ion_intensities": [],
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

            # If line start with "L" and the second field does not start with "decoy_", collect it
            elif line.startswith("L") and not line.split("\t")[1].startswith("decoy_"):
                uniprot_data["uniprot_ids"].append(line.split("\t")[1])

    return uniprot_data


def convert_ids(uniprot_data: dict) -> list[str]:
    """
    This function is responsible for converting a list of uniprot IDs to Gene IDs
    """
    # Create a new BioDBNet object
    biodbnet = BioDBNet()
    # Create an empty list to store the entrez IDs
    entrez_ids: list[str] = []
    uniprot_data["entrez_id"] = []

    pprint(uniprot_data)
    print(f"{len(uniprot_data['uniprot_ids'])} uniprot IDs")
    print(f"{len(uniprot_data['ion_intensities'])} ion intensities")
    exit(2)
    # Loop through the uniprot IDs
    for uniprot_id in uniprot_data:
        # Try to get the entrez ID using BioDBNet.db2db()
        try:
            entrez_id = biodbnet.db2db(
                input_db="UniProt Accession",
                output_db="Gene ID",
                input_values=uniprot_id,
            )
            entrez_id = entrez_id["Gene ID"][0]
        # If there is an error, print it and continue
        except Exception as e:
            print(e)
            continue
        # Add the entrez ID to the list
        entrez_ids.append(entrez_id)
    # Return the list of entrez IDs
    return entrez_ids


def find_sqt_file(input_dir: str) -> set[str]:
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


def main():
    """
    This is the main function for the script
    It is responsible for calling r
    """
    # Define the input directory
    input_dir: str = os.path.join("/Users/joshl/Downloads/b_cell_results")

    # Get a list of sqt files
    sqt_files: set[str] = find_sqt_file(input_dir)

    # Loop through the sqt files and save them to a new list

    read_uniprot_ids(
        "/Users/joshl/Downloads/b_cell_results/PXD026140/comet.20141016_RP-4H_15E5Bcells_120min_top2DD_HCD_01.target.sqt"
    )

    # for sqt_file in sqt_files:
    #     # Get the uniprot IDs from the sqt file
    #     uniprot_ids: dict = read_uniprot_ids(sqt_file)
    #     # Convert uniprot IDs to a new list of entrez IDs
    #     entrez_ids: list[str] = convert_ids(uniprot_ids)


if __name__ == "__main__":
    main()
