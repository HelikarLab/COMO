from bioservices import BioDBNet
import csv
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import os
from pathlib import Path
import subprocess


class RAWtoMZML:
    def __init__(
            self,
            raw_file_input: list[Path],
            mzml_output_dir: Path,
            core_count: int
    ):
        """
        This class is responsible for converting RAW files to mzML format
        """
        self._raw_files: list[Path] = raw_file_input
        self._mzml_output_dir: Path = mzml_output_dir
        self._core_count: int = core_count
        
        # These items are used to track the progress of the conversion
        self._conversion_counter: Synchronized = multiprocessing.Value("i", 1)
        self._finished_counter: Synchronized = multiprocessing.Value("i", 1)
        
        # Create a manager so each process can append data to variables
        # From: https://stackoverflow.com/questions/67974054
        self._mzml_file_paths: list[Path] = multiprocessing.Manager().list()

        # ----- Function Calls -----
        self.raw_to_mzml_wrapper()  # Convert from raw files to mzML

    def raw_to_mzml_wrapper(self):
        """
        This is a wrapper function to multiprocess converting raw files to mzML using ThermoRawFileParser
        """
        file_chunks: list[Path] = np.array_split(self._raw_files, self._core_count)

        jobs: list[multiprocessing.Process] = []

        for i, files in enumerate(file_chunks):
            # Parenthesis + comma needed to make tuple in "args"
            job = multiprocessing.Process(target=self.raw_to_mzml, args=(files,))
            jobs.append(job)

        [job.start() for job in jobs]  # Start jobs
        [job.join() for job in jobs]  # Wait for jobs to finish
        [job.terminate() for job in jobs]  # Terminate jobs
        
    def raw_to_mzml(
        self,
        files: list[Path],
    ):
        """
        This function is responsible or converting the list of raw files to mzML format
        """
        for file in files:

            file_name = f"{file.stem}.mzml"
            save_path: Path = Path(self._mzml_output_dir, file_name)
            self._mzml_file_paths.append(save_path)

            self._conversion_counter.acquire()
            print(f"Starting raw -> mzML conversion: {self._conversion_counter.value} / {len(self._raw_files)} - {file}")  # fmt: skip
            self._conversion_counter.value += 1
            self._conversion_counter.release()

            process = subprocess.run(
                [
                    "thermorawfileparser",
                    "--input",
                    str(file),
                    "--output_file",
                    save_path,
                ],
                stdout=subprocess.PIPE,
            )

            self._finished_counter.acquire()
            print(f"\tFinished raw -> mzML conversion: {self._finished_counter.value} / {len(self._raw_files)} - {file_name}")  # fmt: skip
            self._finished_counter.value += 1
            self._finished_counter.release()

    @property
    def mzml_file_paths(self) -> list[Path]:
        return self._mzml_file_paths
        
        
class MZMLtoSQT:
    def __init__(
        self,
        mzml_file_paths: list[Path],
        fasta_database: Path,
        sqt_output_dir: Path,
        core_count: int,
    ):
        """
        This file is responsible for calling the crux-toolkit utilities to process raw proteomics data
        The tools here are called through the command line

        The following steps must be performed:
        1. Convert RAW files to mzML using thermorawfileparser, saving these to the default mzML output directory
        2. Analyze the mzML files using Crux Comet
        3. Save the SQT files to the default SQT output directory
        """
        # These items are passed into the class
        self._mzml_files: list[Path] = mzml_file_paths
        self._fasta_database: Path = fasta_database
        self._sqt_output_dir: Path = sqt_output_dir
        self._core_count: int = core_count
        
        # Create a manager so each process can append data to variables
        # From: https://stackoverflow.com/questions/67974054
        self._sqt_file_paths: list[Path] = multiprocessing.Manager().list()

        # ----- Function Calls -----
        self.mzml_to_sqt()          # Analyze mzML files, creating SQT files
        
    def mzml_to_sqt(self):
        """
        This function analyzes the converted mzML files and creates SQT files
        This function does not use multiprocessing, as Crux Comet incorporates its own multiprocessing
        """
        for i, file_path in enumerate(self._mzml_files):
            print(f"Creating SQT: {i + 1} / {len(self._mzml_files)} - {file_path.name}")  # fmt: skip
            
            # Call subprocess on command
            subprocess.run(
                [
                    "crux",
                    "comet",
                    "--output_sqtfile", "1",
                    "--output_txtfile", "1",
                    "--output_pepxmlfile", "1",
                    f"--output-dir", self._sqt_output_dir,
                    "--overwrite", "T",
                    "--decoy_search", "1",
                    "--num_threads",
                    str(self._core_count),
                    file_path,
                    self._fasta_database
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            
            # Replace all "comet.*" in output directory with the name of the file being processed
            for file in os.listdir(str(self._sqt_output_dir)):
                if str(file).startswith("comet."):
                    # Get the old file path
                    old_file_path: Path = Path(self._sqt_output_dir, file)
                    
                    # Determine the new file name and path
                    new_file_name = str(file).replace("comet", file_path.stem)
                    new_file_path: Path = Path(self._sqt_output_dir, new_file_name)
                    
                    # Rename the file
                    os.rename(old_file_path, new_file_path)
                    if new_file_path.suffix == ".sqt":
                        self._sqt_file_paths.append(new_file_path)

        print("SQT creation finished")
    
    @property
    def sqt_file_paths(self) -> list[Path]:
        return self._sqt_file_paths


class SQTtoCSV:
    def __init__(self, sqt_input_files: list[Path], output_file: Path):
        """
        This class is meant to convert UniProt IDs to Entrez IDs using BioDBNet
        """
        self._sqt_input_files: list[Path] = sqt_input_files
        self._output_csv: Path = output_file
        
        # Create an empty dictionary to store the UniProt IDs and ion intensities
        # uniprot_ids: dict[str, list[list[str]]
        # ion_intensities: dict[str, list[str]]
        self._uniprot_data: dict = {
            "ion_intensities": [],
            "uniprot_ids": []
        }
        
        self.collect_uniprot_ids_and_ion_intensity()
        self.convert_ids()
        self.write_data()
    
    def _uniprot_from_fasta_header(self, fasta_header: str, separator: str = "|") -> str:
        """
        This function is responsible for collecting the first-index field from a pipe-separated string
        """
        return fasta_header.split(separator)[1]
    
    def collect_uniprot_ids_and_ion_intensity(self):
        """
        This function is responsible for collecting the UniProt IDs from the input sqt files

        Documentation: https://crux.ms/file-formats/sqt-format.html

        We must perform the following
        1. Skip every line starting with an "H"
        2. Skip every line starting with an "M"
        3. Collect every line starting with an "S". This contains the ion intensity
            - Collect the next line that starts with an "L". This contains the UniProt ID
        4. Repeat steps 2 and 3 until the end of the file
        """
        
        for file in self._sqt_input_files:
            
            with open(file, "r") as i_stream:
                for i, line in enumerate(i_stream):
                    
                    # If the line starts with an "S", collect it
                    if line.startswith("S"):
                        # Collect the 8th tab-separated field, it is the ion intensity (index of 7)
                        # Save this in the uniprot_ids dictionary
                        
                        try:
                            ion_intensity = line.split("\t")[7]
                            self._uniprot_data["ion_intensities"].append(ion_intensity)
                        except IndexError:
                            print(f"Error on line {i} in file {file} - Crux collect_uniprot_ids_and_ion_intensity")
                        
                        # Find the next "Locus" line that does not start with "decoy_"
                        locus_found = False
                        while not locus_found:
                            try:
                                line = next(i_stream)
                            except StopIteration:  # We have reached the end of the file
                                break
                            
                            if line[0] == "L":
                                # Skip line if it starts with "decoy_"
                                fasta_header = line.split("\t")[1]
                                if fasta_header.startswith("decoy_"):
                                    continue
                                else:
                                    # Get the uniprot ID from the header
                                    uniprot_id = self._uniprot_from_fasta_header(fasta_header)
                                    self._uniprot_data["uniprot_ids"].append(uniprot_id)
                                    locus_found = True
    
    def convert_ids(self):
        """
        This function converts a list of uniprot IDs to Gene IDs
        """
        biodbnet = BioDBNet()
        self._uniprot_data["gene_ids"] = []
        # Create a progress bar so we know how long this takes
        max_iteration: int = len(self._uniprot_data["uniprot_ids"])
        # for i in tqdm.tqdm(range(0, max_iteration, step)):
        for i in range(0, max_iteration, 500):
            print(f"\rCollecting Gene IDs - {i} to {i + 500} of {max_iteration}", end="", flush=True)
            # Get the next 500 UniProt IDs
            id_lookup = self._uniprot_data["uniprot_ids"][i: i + 500]
            
            # Convert UniProt IDs to Gene IDs
            gene_ids = biodbnet.db2db("UniProt Accession", "Gene ID", input_values=id_lookup)
            
            # Add values to list
            self._uniprot_data["gene_ids"].extend(gene_ids["Gene ID"])
        print("")  # go to next line; we have a carraige return in "for" loop
    
    def write_data(self):
        """
        This function is responsible for writing a dictionary to a csv file

        The dictionary has the following keys:
        1. uniprot_ids
        2. gene_ids
        3. ion_intensities

        Each key has a list of values
        """
        
        with open(self._output_csv, "w") as o_stream:
            writer = csv.writer(o_stream, delimiter=",")
            writer.writerow(["uniprot_id", "gene_id", "ion_intensity"])  # header
            
            # Write data
            for uniprot_id, gene_id, ion_intensity in zip(
                    self._uniprot_data["uniprot_ids"],
                    self._uniprot_data["gene_ids"],
                    self._uniprot_data["ion_intensities"]
            ):
                writer.writerow([uniprot_id, gene_id, ion_intensity])
    
    @property
    def output_csv(self):
        return self._output_csv


if __name__ == "__main__":
    print("Use the main_proteomics.py file")
