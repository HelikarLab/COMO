#!/usr/bin/python3

import os
import tomllib
from datetime import datetime
from dataclasses import dataclass, field
from typing import Literal

@dataclass
class general:
    taxon_id: str
    context_names: list[str] = field(default_factory=list)
    
@dataclass
class rnaseq_preprocess:
    create_counts_matrix: bool
    gene_format: str
    preprocess_mode: Literal["provide-matrix", "create-matrix"]
    
@dataclass
class rna_seq_generation:
    trnaseq_config_file: str
    mrnaseq_config_file: str
    technique: str
    rep_ratio: float
    group_ratio: float
    rep_ratio_h: float
    group_ratio_h: float
    quantile: int
    min_zfpkm: int

@dataclass
class proteomics_analysis:
    proteomics_config_file: str
    rep_ratio: float
    batch_ratio: float
    high_rep_ratio: float
    high_batch_ratio: float
    quantile: int
    
@dataclass
class merge_xomics:
    expression_requirement: int
    requirement_adjust: str
    total_rna_weight: int
    mrna_weight: int
    single_cell_weight: int
    proteomics_weight: int
    
@dataclass
class model_creation:
    low_threshold: int
    high_threshold: int
    output_filetypes: str
    objective_dict: dict
    general_model_file: str
    solver: str
    boundary_reactions_filename: str
    force_reactions_filename: str
    exclude_reactions_filename: str
    recon_algorithms: list[str] = field(default_factory=list[str])

@dataclass
class disease_analysis:
    data_source: str
    disease_names: list[str] = field(default_factory=list[str])

@dataclass
class drug_repurposing:
    sovler: str
    drug_raw_file: str
    
@dataclass
class about:
    version: str = field(init=False)
    date: str = field(init=False)
    
    def __post_init__(self):
        # Set self.date to todays date in the format of "Month Day, Year"
        self.date = datetime.now().strftime("%B %d, %Y")
        
        # Get "VERSION" from the environment
        self.version = os.environ["VERSION"]


class Configs:
    def __init__(self, projectdir):
        self.rootdir = projectdir
        self.datadir = os.path.join(projectdir, "data")
        self.configdir = os.path.join(projectdir, "data", "config_sheets")
        self.outputdir = os.path.join(projectdir, "output")
        self.pydir = os.path.join(projectdir, "py")
       
        self._toml_data: dict = self._read_from_toml()
        self.general = self._toml_data["general"]
        self.rnaseq_preprocess = self._toml_data["rnaseq_preprocess"]
        self.rna_seq_generation = self._toml_data["rna_seq_generation"]
        self.proteomics_analysis = self._toml_data["proteomics_analysis"]
        self.merge_xomics = self._toml_data["merge_xomics"]
        self.model_creation = self._toml_data["model_creation"]
        self.disease_analysis = self._toml_data["disease_analysis"]
        self.drug_repurposing = self._toml_data["drug_repurposing"]
        self.about = self._toml_data["about"]
        
        self._add_date_to_about()

    def _read_from_toml(self):
        toml_file: str = os.path.join(self.rootdir, "config.toml")
        with open(toml_file, "rb") as i_stream:
            data = tomllib.load(i_stream)
        return data
    
    def _add_date_to_about(self):
        """
        This function will add the current date in the format of "Month Day, Year" to the self.about dictionary
        """
        self.about["date"] = datetime.now().strftime("%B %d, %Y")
    

current_dir = os.getcwd()
directory_list = current_dir.split("/")

# Find the "main" directory
split_index = 1
for directory in directory_list:
    if directory == "main":
        break  # Exit the loop when we find the "main" directory
    split_index += 1  # Otherwise increment the index

# Unpack items in dirlist
# From: https://stackoverflow.com/questions/14826888
work_dir = os.path.join(*directory_list[0:split_index])

# Add leading "/", as it will not exist right now
work_dir = os.path.join("/", work_dir)
configs = Configs(work_dir)
