"""
TODO: Integrate crux percolator into this workflow
"""
from bioservices import BioDBNet
from FileInformation import FileInformation
from FileInformation import clear_print
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import os
import pandas as pd
from pathlib import Path
import subprocess


class RAWtoMZML:
    def __init__(self, file_information: list[FileInformation], core_count: int):
        """
        This class is responsible for converting RAW files to mzML format
        """
        self.file_information: list[FileInformation] = file_information
        self._core_count: int = core_count

        # These items are used to track the progress of the conversion
        self._conversion_counter: Synchronized = multiprocessing.Value("i", 0)

        # ----- Function Calls -----
        self.raw_to_mzml_wrapper()  # Convert from raw files to mzML

    def raw_to_mzml_wrapper(self):
        """
        This is a wrapper function to multiprocess converting raw files to mzML using ThermoRawFileParser
        """
        print("Starting raw -> mzML conversion")
        file_chunks: list[Path] = np.array_split(self.file_information, self._core_count)

        jobs: list[multiprocessing.Process] = []
        for i, information in enumerate(file_chunks):
            # Parenthesis + comma needed to make tuple in "args"
            job = multiprocessing.Process(target=self.raw_to_mzml, args=(information,))
            jobs.append(job)

        [job.start() for job in jobs]  # Start jobs
        [job.join() for job in jobs]  # Wait for jobs to finish
        [job.terminate() for job in jobs]  # Terminate jobs

    def raw_to_mzml(self, file_information: list[FileInformation]):
        """
        This function is responsible or converting the list of raw files to mzML format
        """
        for information in file_information:

            self._conversion_counter.acquire()
            self._conversion_counter.value += 1
            clear_print(f"Starting mzML conversion: {self._conversion_counter.value} / {len(self.file_information)} - {information.raw_file_name}")
            self._conversion_counter.release()
            
            information.mzml_base_path.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                [
                    "thermorawfileparser",
                    f"--input={str(information.raw_file_path)}",
                    f"--output_file={str(information.mzml_file_path)}"
                ],
                stdout=subprocess.PIPE,
            )


class MZMLtoSQT:
    def __init__(
        self,
        file_information: list[FileInformation],
        fasta_database: Path,
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
        self._file_information: list[FileInformation] = file_information
        self._fasta_database: Path = fasta_database
        self._core_count: int = core_count

        # ----- Function Calls -----
        self.mzml_to_sqt()  # Analyze mzML files, creating SQT files

    def mzml_to_sqt(self):
        """
        This function analyzes the converted mzML files and creates SQT files
        This function does not use multiprocessing, as Crux Comet incorporates its own multiprocessing
        """
        for i, file_information in enumerate(self._file_information):
            # Clear the previous line of output. Required if the new line is shorter than the previous line
            clear_print(f"Creating SQT: {i + 1} / {len(self._file_information)} - {file_information.sqt_file_name}")

            # Call subprocess on command
            # Only create the SQT file
            subprocess.run(
                [
                    "crux",
                    "comet",
                    "--output_sqtfile",
                    "1",
                    "--output-dir",
                    file_information.sqt_base_path,
                    "--overwrite",
                    "T",
                    "--decoy_search",
                    "1",
                    "--num_threads",
                    str(self._core_count),
                    file_information.mzml_file_path,
                    self._fasta_database,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # Replace all "comet.*" in output directory with the name of the file being processed
            index = 0
            comet_files = [
                str(file)
                for file in os.listdir(file_information.sqt_base_path)
                if str(file).startswith("comet.")
            ]
            current_files = os.listdir(file_information.sqt_base_path)
            for file_name in comet_files:
                file_extension = Path(file_name).suffix

                # Determine the old file path
                old_file_path: Path = Path(file_information.sqt_base_path, file_name)

                # Determine the new file path
                new_file_name = file_name.replace("comet", file_information.base_name)
                new_file_path: Path = Path(
                    file_information.sqt_base_path, new_file_name
                )

                # Rename the file
                os.rename(old_file_path, new_file_path)

        clear_print("SQT creation finished")
        print("")


class SQTtoCSV:
    def __init__(self, file_information: list[FileInformation], core_count: int):
        """
        This class is meant to convert UniProt IDs to Entrez IDs using BioDBNet
        """
        self._file_information: list[FileInformation] = file_information
        self._core_count: int = min(core_count, 4)  # Maximum of 4 cores
        
        # Multiprocessing items
        # File counters: Store the files processed and the total number to process
        # Gene counters: Store the progress of overall genes converted
        self._current_file_counter: Synchronized = multiprocessing.Value("i", 0)
        self._max_file_count: Synchronized = multiprocessing.Value("i", len(self._file_information))
        self._current_gene_counter: Synchronized = multiprocessing.Value("i", 0)
        self._max_gene_count: Synchronized = multiprocessing.Value("i", 0)
        
        self.collect_uniprot_ids_and_ion_intensity()
        self._convert_ids_wrapper()
        self._merge_dataframes()
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
        for i, file_information in enumerate(self._file_information):
            # Create a dictionary with strings as the keys and lists as the values
            # uniprot_id will be a list of strings
            # ion_intensity will be a list of floats
            ion_intensities: list[float] = []
            fasta_headers: list[list[str]] = []

            average_intensities_dict: dict = {
                "uniprot": [],
                file_information.prefix: [],
            }
            # average_intensities: pd.DataFrame = pd.DataFrame(columns=["uniprot", replicate_name])  # fmt: skip

            with open(file_information.sqt_file_path, "r") as i_stream:
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
                    average_intensities_dict[file_information.prefix].append(
                        current_intensity
                    )
            average_intensities_df: pd.DataFrame = pd.DataFrame(average_intensities_dict)  # fmt: skip
            average_intensities_df = average_intensities_df.groupby("uniprot", as_index=False).mean()  # fmt: skip

            # Merge the average intensities to the dataframe
            file_information.intensity_df = pd.merge(
                file_information.intensity_df,
                average_intensities_df,
                on="uniprot",
                how="outer",
            )

            # self._intensities = self._intensities.merge(average_intensities_df, on="uniprot", how="outer")  # fmt: skip
            # self._intensities.join(average_intensities_df)
            # self._intensities = pd.concat([self._intensities, average_intensities_df], axis=0)  # fmt: skip

    def _convert_ids_wrapper(self):
        """
        This function is a multiprocessing wrapper around the convert_ids function
        """
        print("Starting UniProt -> Gene Symbol conversion")
        
        # Split file information to chunks
        information_chunks: list[FileInformation] = np.array_split(self._file_information, self._core_count)
        
        # Create a list of jobs
        jobs: list[multiprocessing.Process] = []
        
        # Set the maximum number of genes to be processed
        self._max_gene_count.value = sum(
                [
                    len(frame.intensity_df)
                    for frame in self._file_information
                ]
            )
        
        for i, information_list in enumerate(information_chunks):
            job = multiprocessing.Process(
                target=self.convert_ids_multiprocess,
                args=(information_list,),
            )
            jobs.append(job)
        
        [job.start() for job in jobs]       # Start jobs
        [job.join() for job in jobs]        # Wait for all jobs to finish
        [job.terminate() for job in jobs]   # Terminate jobs
        
        # A new line is required to ensure we are on the next line after multiprocessing processing
        print("")

    def convert_ids_multiprocess(self, file_information: list[FileInformation]):
        """
        This function is responsible for converting UniProt IDs to Gene Symbols
        """
        # Create a BioDBNet object
        # The objet must be created within each process, as BioDBNet cannot be pickled
        biodbnet: BioDBNet = BioDBNet(verbose=False)
        
        # Iterate over each file information
        for information in file_information:
            # Create appropriate chunk sizes
            chunk_size: int = 500
            num_chunks: int = np.ceil(len(information.intensity_df) / chunk_size)
            chunk_data: pd.DataFrame = np.array_split(information.intensity_df, num_chunks)
            
            lower_iteration: int = 0
            upper_iteration: int = 0
            
            for chunk in chunk_data:
                upper_iteration += len(chunk)
                
                self._current_gene_counter.acquire
                self._current_file_counter.acquire
                
                clear_print(f"Converting to Gene Symbols: {self._current_gene_counter.value + 1:,} of {self._max_gene_count.value:,} (processing file {self._current_file_counter.value + 1} / {self._max_file_count.value})")
                self._current_gene_counter.release
                self._current_file_counter.release
                
                input_data = chunk["uniprot"].values.tolist()
                
                # Convert UniProt IDs to Gene Symbols
                gene_symbols: pd.DataFrame = biodbnet.db2db("UniProt Accession", "Gene Symbol", input_values=input_data)
                
                # Wrangle gene symbols into a format able to be updated with file_intensities
                gene_symbols["uniprot"] = gene_symbols.index
                gene_symbols.rename(columns={"Gene Symbol": "symbol"}, inplace=True)
                
                # Reindexing is the only way to ensure gene_ids are placed at the correct spot
                gene_symbols.reset_index(inplace=True, drop=True)
                gene_symbols.index = range(lower_iteration, upper_iteration)
                
                # Update the file_information.intensity_df with the new gene symbols
                information.intensity_df.update(gene_symbols)
                
                # Update the iterations
                lower_iteration += len(chunk)
                self._current_gene_counter.acquire
                self._current_gene_counter.value += len(chunk)
                self._current_gene_counter.release
                
            # Update the file counter to show the current file is done
            self._current_file_counter.acquire
            self._current_file_counter.value += 1
            self._current_file_counter.release
            
            # Show a final print-statement with all updates
            clear_print(f"Converting to Gene Symbols: {self._current_gene_counter.value:,} of {self._max_gene_count.value:,} (processing file {self._current_file_counter.value} / {self._max_file_count.value})")

    def _merge_dataframes(self):
        """
        This function is responsible for merging dataframes that have the same Gene Symbol
        The resulting UniProt IDs will be combined into a single cell, separated by a semicolon ";"
        The intensity values will be averaged across all UniProt IDs
        """
        
        for file_information in self._file_information:
            data_frame: pd.DataFrame = file_information.intensity_df
            
            # Group each dataframe by the symbol, replacing uniprot column with a semicolon separated list of uniprot IDs
            data_frame = data_frame.groupby("symbol").apply(lambda dataframe: self._handle_merge(dataframe))
            
            # Reset the index to ensure the dataframe is in the correct order
            data_frame.reset_index(inplace=True, drop=True)
            
            # Update the dataframe with the new merged data
            file_information.intensity_df = data_frame
            
    def _handle_merge(self, data_frame: pd.DataFrame):
        
        # Combine all uniprot values in the incoming dataframe into a single string
        # This only combines uniprot values that have the same symbol, as a result of the .apply() function called previously
        uniprot_column: str = ";".join(data_frame["uniprot"].tolist())
        
        # Average the intensity values
        # From: stackoverflow.com/a/32751412/
        new_dataframe: pd.DataFrame = data_frame.groupby("symbol").mean().reset_index()
        
        # Drop "Unnamed: [NUMBER]" columns
        # The "Unnamed" column is formed as a result of the .groupby() call
        # From: https://stackoverflow.com/a/43983654/
        new_dataframe = new_dataframe.loc[:, ~new_dataframe.columns.str.contains("^Unnamed: [0-9]+$")]
        
        # Add the uniprot values back to the dataframe at index 1
        new_dataframe.insert(loc=1, column="uniprot", value=uniprot_column)
        
        return new_dataframe
        
    def write_data(self):
        """
        This function creates a unique dataframe for each cell type found in the intensity dataframes
            from the self._file_information list
        It merges these intensity dataframes, creating a new column for each dataframe within each cell type

        The final dataframes will have the following headers:
        1. uniprot_ids
        2. gene_ids
        3. ion_intensities

        It then writes the dataframes to separate CSV files, dependent on the cell type
        This function is responsible for writing a dictionary to a csv file

        The CSV will be written to the intensity_csv value within each FileInformation object
        """
        # Create a dictionary containing the cell type as keys and the final dataframe as values
        #   This will be used to write the dataframes to separate CSV files
        master_frames: dict[str, pd.DataFrame] = {}

        # Iterate through each FileInformation object
        for file_information in self._file_information:
            # Create a new dataframe for each cell type
            if file_information.cell_type not in master_frames:
                parent_directory: Path = Path(file_information.intensity_csv).parent
                parent_directory.mkdir(parents=True, exist_ok=True)

                master_frames[file_information.cell_type] = pd.DataFrame(
                    columns=file_information.df_columns
                )

            # Update the master frame for the current cell type
            # The master frame should be matched by the uniprot column
            master_frames[file_information.cell_type] = pd.merge(
                master_frames[file_information.cell_type],
                file_information.intensity_df,
                on=["uniprot", "symbol"],
                how="outer",
            )

        # Once merging is complete, write each cell type to its CSV file
        # TODO: Gene Symbols are being repeated as a result of UniProt isoforms. Not sure how to fix this yet
        for file_information in self._file_information:
            master_frames[file_information.cell_type].to_csv(
                file_information.intensity_csv, index=False, na_rep="0"
            )


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
