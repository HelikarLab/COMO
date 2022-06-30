import csv
import os
import sys
from urllib.parse import urlparse
from pathlib import Path
from io import StringIO

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteomics.proteomics_preprocess import ParseCSVInput



def TestParseCSVInput():
    """
    C
    """
    # Create a test CSV file in a temporary location
    file_data: list[list[str]] = [["url,cell_type,study"], ["ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026140,naiveB,S1"]]
    file_path = StringIO()
    csv.writer(file_path).writerows(file_data)
    for line in file_path:
        print(line)
    

if __name__ == '__main__':

