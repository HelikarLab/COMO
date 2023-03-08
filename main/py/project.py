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
        self.version = os.environ["COMO_VERSION"]


class Configs:
    def __init__(self, projectdir):
        self.rootdir = projectdir
        self.datadir = os.path.join(projectdir, "data")
        self.configdir = os.path.join(projectdir, "data", "config_sheets")
        self.outputdir = os.path.join(projectdir, "output")
        self.pydir = os.path.join(projectdir, "py")
       
        self._toml_data: dict = self._read_from_toml()
        
        self.general: general = self._get_general()
        self.rnaseq_preprocess: rnaseq_preprocess = self._get_rnaseq_preprocess()
        self.rna_seq_generation: rna_seq_generation = self._get_rna_seq_generation()
        self.proteomics_analysis: proteomics_analysis = self._get_proteomics_analysis()
        self.merge_xomics: merge_xomics = self._get_merge_xomics()
        self.model_creation: model_creation = self._get_model_creation()
        self.disease_analysis: disease_analysis = self._get_disease_analysis()
        self.drug_repurposing: drug_repurposing = self._get_drug_repurposing()
        self.about: about = self._get_about()

    def _get_general(self) -> general:
        data: general = general(
            taxon_id=self._toml_data["general"]["taxon_id"],
            context_names=self._toml_data["general"]["context_names"],
        )
        return data

    def _get_rnaseq_preprocess(self) -> rnaseq_preprocess:
        data: rnaseq_preprocess = rnaseq_preprocess(
            create_counts_matrix=self._toml_data["rnaseq_preprocess"]["create_counts_matrix"],
            gene_format=self._toml_data["rnaseq_preprocess"]["gene_format"],
            preprocess_mode=self._toml_data["rnaseq_preprocess"]["preprocess_mode"],
        )
        return data

    def _get_rna_seq_generation(self) -> rna_seq_generation:
        data: rna_seq_generation = rna_seq_generation(
            trnaseq_config_file=self._toml_data["rna_seq_generation"]["trnaseq_config_file"],
            mrnaseq_config_file=self._toml_data["rna_seq_generation"]["mrnaseq_config_file"],
            technique=self._toml_data["rna_seq_generation"]["technique"],
            rep_ratio=self._toml_data["rna_seq_generation"]["rep_ratio"],
            group_ratio=self._toml_data["rna_seq_generation"]["group_ratio"],
            rep_ratio_h=self._toml_data["rna_seq_generation"]["rep_ratio_h"],
            group_ratio_h=self._toml_data["rna_seq_generation"]["group_ratio_h"],
            quantile=self._toml_data["rna_seq_generation"]["quantile"],
            min_zfpkm=self._toml_data["rna_seq_generation"]["min_zfpkm"],
        )
        return data

    def _get_proteomics_analysis(self) -> proteomics_analysis:
        data: proteomics_analysis = proteomics_analysis(
            proteomics_config_file=self._toml_data["proteomics_analysis"]["proteomics_config_file"],
            rep_ratio=self._toml_data["proteomics_analysis"]["rep_ratio"],
            batch_ratio=self._toml_data["proteomics_analysis"]["batch_ratio"],
            high_rep_ratio=self._toml_data["proteomics_analysis"]["high_rep_ratio"],
            high_batch_ratio=self._toml_data["proteomics_analysis"]["high_batch_ratio"],
            quantile=self._toml_data["proteomics_analysis"]["quantile"],
        )
        return data

    def _get_merge_xomics(self) -> merge_xomics:
        data: merge_xomics = merge_xomics(
            expression_requirement=self._toml_data["merge_xomics"]["expression_requirement"],
            requirement_adjust=self._toml_data["merge_xomics"]["requirement_adjust"],
            total_rna_weight=self._toml_data["merge_xomics"]["total_rna_weight"],
            mrna_weight=self._toml_data["merge_xomics"]["mrna_weight"],
            single_cell_weight=self._toml_data["merge_xomics"]["single_cell_weight"],
            proteomics_weight=self._toml_data["merge_xomics"]["proteomics_weight"],
        )
        return data

    def _get_model_creation(self) -> model_creation:
        data: model_creation = model_creation(
            low_threshold=self._toml_data["model_creation"]["low_threshold"],
            high_threshold=self._toml_data["model_creation"]["high_threshold"],
            output_filetypes=self._toml_data["model_creation"]["output_filetypes"],
            objective_dict=self._toml_data["model_creation"]["objective_dict"],
            general_model_file=self._toml_data["model_creation"]["general_model_file"],
            solver=self._toml_data["model_creation"]["solver"],
            boundary_reactions_filename=self._toml_data["model_creation"]["boundary_reactions_filename"],
            force_reactions_filename=self._toml_data["model_creation"]["force_reactions_filename"],
            exclude_reactions_filename=self._toml_data["model_creation"]["exclude_reactions_filename"],
            recon_algorithms=self._toml_data["model_creation"]["recon_algorithms"],
        )
        return data

    def _get_disease_analysis(self) -> disease_analysis:
        data: disease_analysis = disease_analysis(
            data_source=self._toml_data["disease_analysis"]["data_source"],
            disease_names=self._toml_data["disease_analysis"]["disease_names"],
        )
        return data

    def _get_drug_repurposing(self) -> drug_repurposing:
        data: drug_repurposing = drug_repurposing(
            sovler=self._toml_data["drug_repurposing"]["sovler"],
            drug_raw_file=self._toml_data["drug_repurposing"]["drug_raw_file"],
        )
        return data

    def _get_about(self) -> about:
        return about()

    def _read_from_toml(self):
        toml_file: str = os.path.join(self.rootdir, "config.toml")
        with open(toml_file, "rb") as i_stream:
            data = tomllib.load(i_stream)
        return data
    

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
