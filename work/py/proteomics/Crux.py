"""
TODO: Integrate crux percolator into this workflow
"""
import pandas as pd
from bioservices import BioDBNet
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import os
from pathlib import Path
import subprocess


class Utils:
    @staticmethod
    def get_cell_type(input_file: Path | str):
        """
        This function will return the "CellType" portion of a file name
    
        Example:
            Input: naiveB_S1R1_experiment_title.raw
            Output: naiveB
        """
        file_name = Path(input_file)
        cell_type = file_name.stem.split("_")[0]
        return cell_type
    
    @staticmethod
    def get_experiment_number(input_file: Path | str):
        """
        This function will return the "S1R1" portion of a file name
    
        Example:
            Input: naiveB_S1R1_experiment_title.raw
            Output: S1R1
        """
        file_name = Path(input_file)
        experiment = file_name.stem.split("_")[1]
        return experiment
    
    @staticmethod
    def get_replicate_number(input_file: Path | str):
        """
        This function will return the "CellType_S1R1" portion of a file name
    
        Example:
            Input: naiveB_S1R1_experiment_title.raw
            Output: naiveB_S1R1
        """
        file_name = Path(input_file)
        prefix = file_name.stem.split("_")
        prefix = f"{prefix[0]}_{prefix[1]}"
        return prefix


class RAWtoMZML:
    def __init__(
        self, raw_file_input: list[Path], mzml_output_dir: Path, core_count: int
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
            print(f"\tFinished raw -> mzML conversion: {self._finished_counter.value} / {len(self._raw_files)}")  # fmt: skip
            self._finished_counter.value += 1
            self._finished_counter.release()

    @property
    def mzml_file_paths(self) -> list[Path]:
        return self._mzml_file_paths


class MZMLtoSQT:
    def __init__(
        self,
        cell_types: list[str],
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
        self._cell_types: list[str] = cell_types
        self._mzml_files: list[Path] = mzml_file_paths
        self._fasta_database: Path = fasta_database
        self._sqt_output_dir: Path = sqt_output_dir
        self._core_count: int = core_count

        # Create a manager so each process can append data to variables
        # From: https://stackoverflow.com/questions/67974054
        self._sqt_file_paths: list[Path] = multiprocessing.Manager().list()

        # ----- Function Calls -----
        self.mzml_to_sqt()  # Analyze mzML files, creating SQT files

    def mzml_to_sqt(self):
        """
        This function analyzes the converted mzML files and creates SQT files
        This function does not use multiprocessing, as Crux Comet incorporates its own multiprocessing
        """
        for i, file_path in enumerate(self._mzml_files):
            print(f"Creating SQT: {i + 1} / {len(self._mzml_files)} - {file_path.name}")  # fmt: skip

            # Call subprocess on command
            # Only create the SQT file
            subprocess.run(
                [
                    "crux",
                    "comet",
                    "--output_sqtfile",
                    "1",
                    f"--output-dir",
                    self._sqt_output_dir,
                    "--overwrite",
                    "T",
                    "--decoy_search",
                    "1",
                    "--num_threads",
                    str(self._core_count),
                    file_path,
                    self._fasta_database,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # Replace all "comet.*" in output directory with the name of the file being processed
            index = 0
            for file in os.listdir(str(self._sqt_output_dir)):
                if str(file).startswith("comet."):
                    replicate_name: str = f"{self._cell_types[index]}_S1R{index + 1}"
                    # Get the old file path
                    old_file_path: Path = Path(self._sqt_output_dir, file)

                    # Determine the new file name and path
                    # Example: naiveB_S1R1_experiment_name.sqt
                    new_file_name = str(file).replace("comet", f"{replicate_name}_{file_path.stem}")  # fmt: skip
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
    def __init__(self, cell_types: list[str], sqt_input_files: list[Path], output_file: Path):
        """
        This class is meant to convert UniProt IDs to Entrez IDs using BioDBNet
        """
        self._cell_types: list[str] = cell_types
        self._sqt_input_files: list[Path] = sqt_input_files
        self._output_csv: Path = output_file
        self._intensities: pd.DataFrame = pd.DataFrame(
            columns=["uniprot", "symbol"]
        )

        self.collect_uniprot_ids_and_ion_intensity()
        self.convert_ids()
        self.write_data()

    def _uniprot_from_fasta_header(
        self, fasta_header: str, separator: str = "|"
    ) -> str:
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
        for i, file in enumerate(sorted(self._sqt_input_files)):
            # Create a dictionary with strings as the keys and lists as the values
            # uniprot_id will be a list of strings
            # ion_intensity will be a list of floats
            replicate_name: str = f"{Utils.get_replicate_number(file)}_intensity"  # naiveB_S1R1_intensity
            ion_intensities: list[float] = []
            fasta_headers: list[list[str]] = []
            
            average_intensities_dict: dict = {"uniprot": [], replicate_name: []}
            # average_intensities: pd.DataFrame = pd.DataFrame(columns=["uniprot", replicate_name])  # fmt: skip

            with open(file, "r") as i_stream:
                """
                We are going to use spectra_line_nums if the list starts with "S"
                Beneath this, we are going to collect every locus ("L") that does not have a "decoy_" in a nested list
                The result will be each spectra value corresponds to a list of UniProt IDs
                """
                for j, line in enumerate(i_stream):

                    # If the line starts with an "S", collect it
                    if line.startswith("S"):
                        # If the length of ion_intensities is not equal to fasta_headers,
                        #   We have added an intensity that did not have valid locus data
                        # (i.e., it only contained "decoy_")
                        if len(ion_intensities) != len(fasta_headers):
                            ion_intensities.pop()
                        
                        intensity: float = float(line.split("\t")[7])
                        ion_intensities.append(intensity)
                        fasta_headers.append([])

                    # Get sequential lines starting with "L" and append the uniprot ID to the dataframe
                    elif line.startswith("L"):
                        fasta_header = line.split("\t")[1]
                        if fasta_header.startswith("decoy_"):
                            continue
                        else:
                            # Append fasta header to most recent list created
                            fasta_headers[-1].append(fasta_header)

            # Append corresponding values in ion_intensities and fasta_headers to the average_intensities list
            # concat_dict: dict = {"uniprot": [], replicate_name: []}
            for j in range(len(ion_intensities)):
                current_intensity = ion_intensities[j]

                for k in range(len(fasta_headers[j])):
                    current_fasta_header = fasta_headers[j][k]

                    # Get the UniProt ID from the fasta header
                    uniprot_id = self._uniprot_from_fasta_header(current_fasta_header)

                    # Create a new row in the dataframe
                    average_intensities_dict["uniprot"].append(uniprot_id)
                    average_intensities_dict[replicate_name].append(current_intensity)

            # Group by uniprot_id and take the mean of the ion_intensities
            average_intensities_df: pd.DataFrame = pd.DataFrame(average_intensities_dict)  # fmt: skip
            average_intensities_df = average_intensities_df.groupby("uniprot", as_index=False).mean()  # fmt: skip

            # Merge the average intensities to the dataframe
            self._intensities = self._intensities.merge(average_intensities_df, on="uniprot", how="outer")  # fmt: skip
            # self._intensities.join(average_intensities_df)
            # self._intensities = pd.concat([self._intensities, average_intensities_df], axis=0)  # fmt: skip

    def convert_ids(self):
        """
        This function converts a list of uniprot IDs to Gene IDs
        """
        biodbnet = BioDBNet()

        """
        len(DATAFRAME) / X = 500
        X = len(DATAFRAME) / 500
        """

        # Get "chunk_sizes" rows at a time
        # This is not perfect, as num_sections rounds up to the nearest integer
        chunk_size: int = 500
        num_sections: int = np.ceil(len(self._intensities) / chunk_size)
        chunk_data: pd.DataFrame = np.array_split(self._intensities, num_sections)
        lower_iteration: int = 0
        upper_iteration: int = 0
        
        for i, chunk in enumerate(chunk_data):
            upper_iteration += len(chunk)
            print(f"\rConverting to Gene Symbols {lower_iteration + 1}:{upper_iteration} of {len(self._intensities)}", end="", flush=True)
            
            input_data = chunk["uniprot"].values.tolist()

            # Convert UniProt IDs to Gene IDs
            gene_ids: pd.DataFrame = biodbnet.db2db("UniProt Accession", "Gene Symbol", input_values=input_data)
            
            # Wrangle gene_ids into a format able to be updated with self._intensities
            gene_ids["uniprot"] = gene_ids.index
            gene_ids.rename(columns={"Gene Symbol": "symbol"}, inplace=True)
            gene_ids.reset_index(inplace=True, drop=True)
            
            # Reindexing was the only way to ensure gene_ids were placed in the correct spot
            gene_ids.index = range(lower_iteration, upper_iteration)
            # Update intensities dataframe
            self._intensities.update(gene_ids)
            
            # Update the iteration to start at
            lower_iteration += len(chunk)
            
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
        self._intensities.to_csv(self._output_csv, index=False, header=True, na_rep="0")

    @property
    def output_csv(self):
        return self._output_csv


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
