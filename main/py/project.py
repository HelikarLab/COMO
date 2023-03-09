import os
import toml
from enum import Enum
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field

from async_bioservices.input_database import InputDatabase
from async_bioservices.taxon_ids import TaxonIDs


@dataclass
class general:
    taxon_id: TaxonIDs
    context_names: list[str] = field(default_factory=list)
    
    
class PreprocessMode(Enum):
    CREATE_MATRIX = "create-matrix"
    PROVIDE_MATRIX = "provide-matrix"

@dataclass
class rnaseq_preprocess:
    create_counts_matrix: bool
    gene_format: InputDatabase
    preprocess_mode: PreprocessMode
    matrix_filename: Path
    
@dataclass
class rna_seq_generation:
    trnaseq_config_filepath: Path
    mrnaseq_config_filepath: Path
    single_cell_config_filepath: Path
    technique: str
    replicate_ratio: float
    group_ratio: float
    high_replicate_ratio: float
    high_group_ratio: float
    quantile: int
    min_zfpkm: int
    min_cpm_count: int | str

@dataclass
class proteomics_analysis:
    proteomics_config_file: Path
    rep_ratio: float
    batch_ratio: float
    high_rep_ratio: float
    high_batch_ratio: float
    quantile: int
    
@dataclass
class merge_xomics:
    custom_requirement_filepath: Path
    microarray_config_filepath: Path
    expression_requirement: int
    requirement_adjust: str
    total_rna_weight: int
    mrna_weight: int
    single_cell_weight: int
    proteomics_weight: int
    no_high_confidence: bool
    no_NA: bool
    merge_distribution: bool
    keep_gene_score: bool
    
@dataclass
class model_creation:
    low_threshold: int
    high_threshold: int
    output_filetypes: str
    objective_dict: dict
    general_model_file: Path
    solver: str
    boundary_reactions_filepath: Path
    force_reactions_filepath: Path
    exclude_reactions_filepath: Path
    recon_algorithms: list[str] = field(default_factory=list[str])

@dataclass
class disease_analysis:
    data_source: str
    disease_names: list[str] = field(default_factory=list[str])

@dataclass
class drug_repurposing:
    sovler: str
    drug_raw_filepath: Path
    
@dataclass
class about:
    version: str = field(init=False)
    date: str = field(init=False)
    
    def __post_init__(self):
        # Set self.date to todays date in the format of "Month Day, Year"
        self.date = datetime.now().strftime("%B %d, %Y")
        
        # EXAMPLE: get "v1.0.0-BRANCH" from "refs/tags/v1.0.0-BRANCH"
        try:
            self.version = os.environ["COMO_VERSION"].split("/")[-1]
        except KeyError:
            self.version = "UNKNOWN"


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

    def _read_from_toml(self):
        toml_file: str = os.path.join(self.rootdir, "config.toml")
        data: dict = toml.loads(open(toml_file).read())
        return data

    def _get_general(self) -> general:
        taxon_id = self._toml_data["general"]["taxon_id"]
        context_names: list[str] = self._toml_data["general"]["context_names"]
        
        if isinstance(taxon_id, str):
            if taxon_id.lower() in ["human", "homo sapiens"]:
                taxon_id = TaxonIDs.HOMO_SAPIENS
            elif taxon_id.lower() in ["mouse", "mus musculus"]:
                taxon_id = TaxonIDs.MUS_MUSCULUS
            else:
                raise ValueError("The taxon_id setting under 'general' must be either an integer, 'human' or 'mouse'.\nPlease edit `config.toml`")
        
        if not isinstance(context_names, list):
            raise ValueError("The context_names setting under 'general' must be a list, such as `['type1', 'type2']`.\nPlease edit `config.toml`")
        
        data: general = general(
            taxon_id=taxon_id,
            context_names=context_names,
        )
        return data

    def _get_rnaseq_preprocess(self) -> rnaseq_preprocess:
        create_counts_matrix = self._toml_data["rnaseq_preprocess"]["create_counts_matrix"]
        gene_format: str = self._toml_data["rnaseq_preprocess"]["gene_format"].lower()
        preprocess_mode: str = self._toml_data["rnaseq_preprocess"]["preprocess_mode"].lower()
        matrix_filename = self._toml_data["rnaseq_preprocess"]["matrix_filename"]
        
        if not isinstance(create_counts_matrix, bool):
            raise ValueError("The create_counts_matrix setting under 'rnaseq_preprocess' must be either 'true' or 'false'.\nPlease edit `config.toml`")
        
        if gene_format not in ["entrez", "ensembl", "symbol"]:
            raise ValueError("The gene_format setting under 'rnaseq_preprocess' must be either 'Entrez', 'Ensembl', or 'Symbol'.\nPlease edit `config.toml`")
        else:
            if gene_format in ["ensembl", "ensemble", "ensg", "ensmusg", "ensembl id", "ensembl gene id"]:
                gene_format_database: InputDatabase = InputDatabase.ENSEMBL_GENE_ID
            elif gene_format in ["hgnc symbol", "hugo", "hugo symbol", "symbol", "hgnc", "gene symbol"]:
                gene_format_database: InputDatabase = InputDatabase.GENE_SYMBOL
            elif gene_format in ["entrez", "entres", "entrez id", "entrez number" "gene id"]:
                gene_format_database: InputDatabase = InputDatabase.GENE_ID
        
        if preprocess_mode not in ["provide-matrix", "create-matrix"]:
            raise ValueError("The preprocess_mode setting under 'rnaseq_preprocess' must be either 'provide-matrix' or 'create-matrix'.\nPlease edit `config.toml`")
        elif preprocess_mode == "provide-matrix" and matrix_filename == "":
            raise ValueError("The matrix_filename setting under 'rnaseq_preprocess' must be set if the preprocess_mode is set to 'create-matrix'.\nPlease edit `config.toml`")
        else:
            if preprocess_mode == "provide-matrix":
                preprocess_mode: PreprocessMode = PreprocessMode.PROVIDE_MATRIX
            elif preprocess_mode == "create-matrix":
                preprocess_mode: PreprocessMode = PreprocessMode.CREATE_MATRIX
        
        data: rnaseq_preprocess = rnaseq_preprocess(
            create_counts_matrix=create_counts_matrix,
            gene_format=gene_format_database,
            preprocess_mode=preprocess_mode,
            matrix_filename=matrix_filename,
        )
        return data

    def _get_rna_seq_generation(self) -> rna_seq_generation:
        scrnaseq_file = self._toml_data["rna_seq_generation"]["single_cell_config_file"]
        if scrnaseq_file == "None":
            scrnaseq_config_filepath: Path = Path()
        else:
            scrnaseq_config_filepath: Path = Path(
                self.configdir,
                self._toml_data["rna_seq_generation"]["single_cell_config_file"]
            )

        trnaseq_config_filepath: Path = Path(self.configdir, self._toml_data["rna_seq_generation"]["trnaseq_config_file"])
        mrnaseq_config_filepath: Path = Path(self.configdir, self._toml_data["rna_seq_generation"]["mrnaseq_config_file"])
        technique: str = self._toml_data["rna_seq_generation"]["technique"]
        rep_ratio: float = self._toml_data["rna_seq_generation"]["replicate_ratio"]
        group_ratio: float = self._toml_data["rna_seq_generation"]["group_ratio"]
        rep_ratio_h: float = self._toml_data["rna_seq_generation"]["high_replicate_ratio"]
        group_ratio_h: float = self._toml_data["rna_seq_generation"]["high_group_ratio"]
        quantile: int = self._toml_data["rna_seq_generation"]["quantile"]
        min_zfpkm: int = self._toml_data["rna_seq_generation"]["min_zfpkm"]
        min_cpm_count = self._toml_data["rna_seq_generation"]["min_cpm_count"]
        
        
        
        if technique.lower() not in ["quantile", "zfpkm", "cpm"]:
            raise ValueError("The technique setting under 'rna_seq_generation' must be either 'quantile', 'zfpkm', or 'cpm'.\nPlease edit `config.toml`")
        
        if min_cpm_count == 0:
            min_cpm_count = "default"

        if int(quantile) > 100 or int(quantile) < 1:
            raise ValueError("Quantile must be between 1 and 100.\nPlease edit `config.toml`")

        data: rna_seq_generation = rna_seq_generation(
            trnaseq_config_filepath=trnaseq_config_filepath,
            mrnaseq_config_filepath=mrnaseq_config_filepath,
            single_cell_config_filepath=scrnaseq_config_filepath,
            technique=technique,
            replicate_ratio=rep_ratio,
            group_ratio=group_ratio,
            high_replicate_ratio=rep_ratio_h,
            high_group_ratio=group_ratio_h,
            quantile=quantile,
            min_zfpkm=min_zfpkm,
            min_cpm_count=min_cpm_count,
        )
        return data

    def _get_proteomics_analysis(self) -> proteomics_analysis:
        proteomics_file = self._toml_data["proteomics_analysis"]["proteomics_config_file"]
        if proteomics_file == "None":
            proteomics_config_file: Path = Path()
        else:
            proteomics_config_file: Path = Path(self.configdir, self._toml_data["proteomics_analysis"]["proteomics_config_file"])
        
        rep_ratio = self._toml_data["proteomics_analysis"]["rep_ratio"]
        batch_ratio = self._toml_data["proteomics_analysis"]["batch_ratio"]
        high_rep_ratio = self._toml_data["proteomics_analysis"]["high_rep_ratio"]
        high_batch_ratio = self._toml_data["proteomics_analysis"]["high_batch_ratio"]
        quantile = self._toml_data["proteomics_analysis"]["quantile"]
        
        data: proteomics_analysis = proteomics_analysis(
            proteomics_config_file=proteomics_config_file,
            rep_ratio=rep_ratio,
            batch_ratio=batch_ratio,
            high_rep_ratio=high_rep_ratio,
            high_batch_ratio=high_batch_ratio,
            quantile=quantile,
        )
        return data

    def _get_merge_xomics(self) -> merge_xomics:
        microarray_file = self._toml_data["merge_xomics"]["microarray_config_file"]
        custom_requirement_file = self._toml_data["merge_xomics"]["custom_requirement_file"]
        if microarray_file == "None":
            microarray_config_filepath: Path = Path()
        else:
            microarray_config_filepath: Path = Path(
                self.configdir,
                self._toml_data["merge_xomics"]["microarray_config_file"]
            )
            
        if custom_requirement_file == "None":
            custom_requirement_filepath: Path = Path()
        else:
            custom_requirement_filepath: Path = Path(
                self.configdir,
                self._toml_data["merge_xomics"]["custom_requirement_file"]
            )

        expression_requirement: int = self._toml_data["merge_xomics"]["expression_requirement"]
        requirement_adjust = self._toml_data["merge_xomics"]["requirement_adjust"]
        total_rna_weight = self._toml_data["merge_xomics"]["total_rna_weight"]
        mrna_weight = self._toml_data["merge_xomics"]["mrna_weight"]
        single_cell_weight = self._toml_data["merge_xomics"]["single_cell_weight"]
        proteomics_weight = self._toml_data["merge_xomics"]["proteomics_weight"]
        no_high_confidence = self._toml_data["merge_xomics"]["no_high_confidence"]
        no_NA = self._toml_data["merge_xomics"]["no_NA"]
        merge_distribution = self._toml_data["merge_xomics"]["merge_distribution"]
        keep_gene_score = self._toml_data["merge_xomics"]["keep_gene_score"]
        
        if expression_requirement < 1:
            raise ValueError("The expression requirement must be at least 1.\nPlease edit `config.toml`")
        
        if requirement_adjust not in ["progressive", "regressive", "flat", "custom"]:
            raise ValueError("The requirement adjust setting under 'merge_xomics' must be either 'progressive', 'regressive', 'flat', or 'custom'.\nPlease edit `config.toml`")
        
        data: merge_xomics = merge_xomics(
            custom_requirement_filepath=custom_requirement_filepath,
            microarray_config_filepath=microarray_config_filepath,
            expression_requirement=expression_requirement,
            requirement_adjust=requirement_adjust,
            total_rna_weight=total_rna_weight,
            mrna_weight=mrna_weight,
            single_cell_weight=single_cell_weight,
            proteomics_weight=proteomics_weight,
            no_high_confidence=no_high_confidence,
            no_NA=no_NA,
            merge_distribution=merge_distribution,
            keep_gene_score=keep_gene_score
        )
        return data

    def _get_model_creation(self) -> model_creation:
        low_threshold = self._toml_data["model_creation"]["low_threshold"]
        high_threshold = self._toml_data["model_creation"]["high_threshold"]
        output_filetypes = self._toml_data["model_creation"]["output_filetypes"]
        objective_dict = self._toml_data["model_creation"]["objective_dict"]
        general_model_file = self._toml_data["model_creation"]["general_model_file"]
        solver = self._toml_data["model_creation"]["solver"]
        boundary_reactions_filepath = self._toml_data["model_creation"]["boundary_reactions_filename"]
        force_reactions_filepath = self._toml_data["model_creation"]["force_reactions_filename"]
        exclude_reactions_filepath = self._toml_data["model_creation"]["exclude_reactions_filename"]
        recon_algorithms = self._toml_data["model_creation"]["recon_algorithms"]
        
        data: model_creation = model_creation(
            low_threshold=low_threshold,
            high_threshold=high_threshold,
            output_filetypes=output_filetypes,
            objective_dict=objective_dict,
            general_model_file=general_model_file,
            solver=solver,
            boundary_reactions_filepath=boundary_reactions_filepath,
            force_reactions_filepath=force_reactions_filepath,
            exclude_reactions_filepath=exclude_reactions_filepath,
            recon_algorithms=recon_algorithms,
        )
        return data

    def _get_disease_analysis(self) -> disease_analysis:
        data_source = self._toml_data["disease_analysis"]["data_source"]
        disease_names = self._toml_data["disease_analysis"]["disease_names"]
        
        data: disease_analysis = disease_analysis(
            data_source=data_source,
            disease_names=disease_names,
        )
        return data

    def _get_drug_repurposing(self) -> drug_repurposing:
        sovler = self._toml_data["drug_repurposing"]["sovler"]
        drug_raw_filepath = self._toml_data["drug_repurposing"]["drug_raw_file"]
        
        data: drug_repurposing = drug_repurposing(
            sovler=sovler,
            drug_raw_filepath=drug_raw_filepath,
        )
        return data

    def _get_about(self) -> about:
        return about()
    

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
