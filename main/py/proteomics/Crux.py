"""
TODO: Integrate crux percolator into this workflow
"""
import asyncio
from bioservices import BioDBNet
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import subprocess
import tqdm

from .FileInformation import FileInformation
from .FileInformation import clear_print


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
        print("")  # Get a new line to print on
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
                    "crux", "comet",
                    "--output_sqtfile", "1",
                    "--output-dir", file_information.sqt_base_path,
                    "--overwrite", "T",
                    "--decoy_search", "1",
                    "--num_threads", str(self._core_count),
                    file_information.mzml_file_path,  # Input mzML file
                    self._fasta_database,  # Database to search
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # Replace all "comet.*" in output directory with the name of the file being processed
            comet_files = [
                str(file)
                for file in os.listdir(file_information.sqt_base_path)
                if str(file).startswith("comet.")
            ]
            for file_name in comet_files:

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
        self._biodbnet: BioDBNet = BioDBNet(verbose=False)
        self._file_information: list[FileInformation] = file_information
        self._core_count: int = min(core_count, 4)  # Maximum of 4 cores
        
        # Merged frames contains all dataframes
        # Split frames contains the S1/S2/etc dataframes, extracted from Merged frames
        self._merged_frames: dict[str, pd.DataFrame] = {}
        self._split_frames: dict[str, [pd.DataFrame]] = {}

        # Max of 15 asynchronous tasks at once
        # From: https://stackoverflow.com/a/48256949
        # self._semaphore = asyncio.Semaphore(50)
        
        self.collect_uniprot_ids_and_ion_intensity()
        asyncio.run(self._convert_uniprot_wrapper())
        self.create_merged_frame()
        self.new_write_data()

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

            # Use dictionary comprehension to create data dictionary
            average_intensities_dict: dict = {
                key: [] for key in list(file_information.intensity_df.columns)
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
                    average_intensities_dict[file_information.batch].append(current_intensity)
            
            # Fill the "symbol" list with NA, as they are not converted yet
            average_intensities_dict["symbol"] = [np.nan] * len(average_intensities_dict["uniprot"])
            
            # Assign the file_information intensity dataframe to the gathered values
            self._file_information[i].intensity_df = pd.DataFrame(average_intensities_dict)
            self._file_information[i].intensity_df = self._file_information[i].intensity_df.groupby("uniprot", as_index=False).mean()

    async def _convert_uniprot_wrapper(self):
        """
        This function is a multiprocessing wrapper around the convert_ids function
        """
        values = [
            self.async_convert_uniprot(self._file_information[i])
            for i in range(len(self._file_information))
        ]
        
        # Create a progress bar of results
        # From: https://stackoverflow.com/a/61041328/
        progress_bar = tqdm.tqdm(desc="Starting UniProt to Gene Symbol conversion... ", total=len(self._file_information))
        for i, result in enumerate(asyncio.as_completed(values)):
            await result  # Get result from asyncio.as_completed
            progress_bar.set_description(f"Working on {i + 1} of {len(self._file_information)}")
            progress_bar.update()
    
    async def async_convert_uniprot(self, file_information: FileInformation):
        
        chunk_size: int = 400
        num_chunks: int = np.ceil(len(file_information.intensity_df) / chunk_size)
        frame_chunks: pd.DataFrame = np.array_split(file_information.intensity_df, num_chunks)
        
        lower_iteration: int = 0
        upper_iteration: int = 0
        
        for chunk in frame_chunks:
            upper_iteration += len(chunk)
            input_values: list[str] = list(chunk["uniprot"])

            # Limit number of asynchronous calls to value defined in self._semaphore
            loop = asyncio.get_event_loop()
            # async with self._semaphore:
            #     gene_symbols: pd.DataFrame = await loop.run_in_executor(None, self._biodbnet.db2db, "UniProt Accession", "Gene Symbol", input_values)
            gene_symbols: pd.DataFrame = await loop.run_in_executor(None, self._biodbnet.db2db, "UniProt Accession", "Gene Symbol", input_values)
            
            # The index is UniProt IDs. Create a new column of these values
            gene_symbols["uniprot"] = gene_symbols.index
            gene_symbols.rename(columns={"Gene Symbol": "symbol"}, inplace=True)
            
            # Create a new "index" column to reset the index of gene_symbols
            gene_symbols["index"] = range(lower_iteration, upper_iteration)
            gene_symbols.set_index("index", inplace=True, drop=True)
            gene_symbols = pd.merge(gene_symbols, chunk[[file_information.batch]], left_index=True, right_index=True)
            
            lower_iteration += len(chunk)
            file_information.intensity_df.update(gene_symbols)
        
    def create_merged_frame(self):
        """
        This function is responsible for merging all dataframes of a specific cell type into a single master frame
        This will allow for aggregating S1R1/S1R2/etc. dataframes into a single S1 dataframe
        """
        for file in self._file_information:
            cell_type = file.cell_type
            if cell_type not in self._merged_frames.keys():
                self._merged_frames[cell_type] = pd.DataFrame(columns=["symbol", "uniprot"])
        
            self._merged_frames[cell_type] = pd.concat([self._merged_frames[cell_type], file.intensity_df])

            # Drop the 'uniprot' column and merge by cell type for each dataframe in master_frame
            self._merged_frames[cell_type].drop(columns=["uniprot"], inplace=True)
            self._merged_frames[cell_type] = self._merged_frames[cell_type].groupby("symbol").mean()
            
            # Create a new column "symbol" that is the index
            self._merged_frames[cell_type]["symbol"] = self._merged_frames[cell_type].index
            
            # Reset the index to ensure the dataframe is in the correct order
            self._merged_frames[cell_type].reset_index(inplace=True, drop=True)
            
            # Replace all nan values with 0
            self._merged_frames[cell_type].fillna(0, inplace=True)
    
    def split_abundance_values(self):
        """
        This function is responsible for splitting abundance values into separate columns based on their S#R# identifier
        
        It will start by finding all S1R*, S2R*, etc. columns in each cell type under self._master_frames
        From here, it will merge these dataframes into a new dataframe under self._split_frames, corresponding to the cell type
        These split frames can then be written by self.write_data()
        
        Example Dataframe:
            Starting
            --------
            symbol,naiveB_S1R1,naiveB_S1R2,naiveB_S2R1
            A     ,100        ,0          ,50
            B     ,200        ,75         ,100
            C     ,150        ,100        ,175
            
            Ending
            ------
            symbol,naiveB_S1,naiveB_S2
            A     ,50       ,50
            B     ,137.5    ,100
            C     ,125      ,175
        """
        # Get a new line to print output on
        print("")
        
        # Copy the dictionary keys from merged_frames into the split_frames
        for key in self._merged_frames.keys():
            if key not in self._split_frames.keys():
                self._split_frames[key] = []
        
        # 'cell_type' is a dictionary key
        for cell_type in self._merged_frames:
            # Must collect the maximum S# value found in the master_frame
            # -------------------
            max_iteration: int = 0
            dataframe = self._merged_frames[cell_type]
            for column in dataframe.columns:
                # Find the {cell_type}_S#R# column using regex
                if re.match(rf"{cell_type}_S\d+R\d+", column):
                    # Find the S# value using regex
                    iteration: int = int(re.search(r"S(\d+)", column).group(1))
                    
                    # Find the maximum S# value
                    if iteration > max_iteration:
                        max_iteration = iteration
            # -------------------
            
            # Aggregate all S# values from 1 to max_iteration
            # The new dataframes will go into the self._split_frames
            for i in range(1, max_iteration + 1):
                # Create a new dataframe to split the S# columns from
                split_frame: pd.DataFrame = dataframe.copy()
                # Get the current S{i} columns in
                abundance_columns: list[str] = [column for column in split_frame.columns if re.match(rf"{cell_type}_S{i}R\d+", column)]
                take_columns: list[str] = ["symbol"] + abundance_columns
                average_intensity_name: str = f"{cell_type}_S{i}"
                
                # Calculate average intensities and assign a new column
                average_intensity_values = split_frame[take_columns].mean(axis=1)
                # split_frame.loc[:, take_columns].mean(axis=1)
                split_frame[average_intensity_name] = average_intensity_values

                # Purge the S#R## column names, they are no longer required
                # We now have a new dataframe with "symbol" and "{cell_type}_S{i}" columns
                # Unpack abundance_columns to create a single list
                split_frame.drop(columns=abundance_columns, inplace=True)
            
                # If the "{cell_type}_S{i}" column is 0, remove it
                split_frame = split_frame[split_frame[average_intensity_name] != 0]
                split_frame.reset_index(inplace=True, drop=True)
                
                # Find duplicate "symbol" values and average across them
                # Duplicate symbols are a result of protein isoforms mapping to the same gene
                split_frame = split_frame.groupby("symbol", as_index=False).mean()
                
                self._split_frames[cell_type].append(split_frame)

    def new_write_data(self):
        """
        This function is responsible for writing the dataframes found in self._merge_frames to the respective cell type file
        """
        # Get a list of CSV file locations
        csv_file_location: dict[str, Path] = {}
        for information in self._file_information:
            if information.cell_type not in csv_file_location:
                csv_file_location[information.cell_type] = information.intensity_csv
                
        # Sort columns of each cell type dataframe
        for key in self._merged_frames.keys():
            # Get the "symbol" column so it can be placed at index 0
            symbol_column = self._merged_frames[key].pop("symbol")
            
            # Sort {cell_type}_S# columns
            col_names: list[str] = list(self._merged_frames[key].columns)
            self._merged_frames[key].reindex(sorted(col_names), axis=1)
            
            # Place the "symbol" column back in at index 0
            self._merged_frames[key].insert(0, "symbol", symbol_column)
            
            # Write the dataframe to its appropriate location
            self._merged_frames[key].to_csv(csv_file_location[key], index=False)
            
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

                master_frames[file_information.cell_type] = pd.DataFrame(columns=file_information.base_columns)  # fmt: skip

            # Update the master frame for the current cell type
            # The master frame should be matched by the uniprot column
            master_frames[file_information.cell_type] = pd.merge(
                master_frames[file_information.cell_type],
                file_information.intensity_df,
                on=["uniprot", "symbol"],
                how="outer",
            )

        # Once merging is complete, write each cell type to its CSV file
        for cell_type in master_frames:
            master_frames[cell_type].replace(np.nan, 0, inplace=True)
            master_frames[cell_type].sort_values(by="symbol", inplace=True, ignore_index=True)
            
            csv_path = FileInformation.intensity_file_path(cell_type=cell_type)
            master_frames[cell_type].to_csv(csv_path, index=False)
            
    @property
    def file_information(self):
        """
        Reassigning the values from the incoming file_manager into a shared-memory variable means
            we must provide an option to return the file_information list
        """
        return self._file_information


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
