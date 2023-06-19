"""
This is a python implementation of the rnaseq.R file
"""
import re
import numpy as np
import pandas as pd
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Union
from dataclasses import dataclass, field
import sklearn


class NormalizationTechnique(Enum):
    CPM = "cpm"
    TPM = "tpm"


class PreparationMethod(Enum):
    TOTAL = "total"
    MRNA = "mrna"


@dataclass
class FilterSettings:
    technique: str
    min_count: int
    min_zfpkm: int
    quantile: float
    batch_ratio: float
    replicate_ratio: float
    batch_ratio_high: float
    replicate_ratio_high: float
    
    def __post_init__(self):
        valid_techniques: list[str] = ["cpm", "zfpkm", "quantile", "umi"]
        self.technique = self.technique.lower()
        if self.technique not in valid_techniques:
            raise ValueError(f"Invalid technique: {self.technique}. Must be one of {valid_techniques}")


@dataclass
class SampleMetricItem:
    study_number: str
    layout: str
    sample_names: str
    count_matrix: pd.DataFrame
    fragment_lengths: int
    entrez_ids: pd.Series
    entrez_ids_high_confidence: pd.Series
    num_samples: int
    gene_size: pd.Series
    z_score: pd.DataFrame
    fpkm_matrix: pd.DataFrame
    tpm_matrix: pd.DataFrame
    cpm_matrix: pd.DataFrame
    
    def __post_init__(self):
        if self.layout not in ["paired-end", "single-end"]:
            raise ValueError(f"Invalid layout: {self.layout}. Must be one of ['paired-end', 'single-end']")


@dataclass
class SampleMetrics:
    study_number: str = None
    layout: list[str] = field(default_factory=list[str])
    sample_names: list[str] = field(default_factory=list[str])
    count_matrix: pd.DataFrame = None
    fragment_lengths: list[int] = field(default_factory=list[int])
    entrez_ids: pd.Series = None
    entrez_ids_high_confidence: pd.Series = None
    num_samples: int = None
    gene_size: pd.Series = None
    z_score: pd.DataFrame = None
    fpkm_matrix: pd.DataFrame = None
    tpm_matrix: pd.DataFrame = None
    cpm_matrix: pd.DataFrame = None
    
    def __getitem__(self, index: Union[int, slice]) -> SampleMetricItem:
        return SampleMetricItem(
            study_number=self.study_number,
            layout=self.layout[index].lower(),
            sample_names=self.sample_names[index],
            count_matrix=self.count_matrix,
            fragment_lengths=self.fragment_lengths[index],
            entrez_ids=self.entrez_ids,
            entrez_ids_high_confidence=self.entrez_ids_high_confidence,
            num_samples=self.num_samples,
            gene_size=self.gene_size,
            z_score=self.z_score,
            fpkm_matrix=self.fpkm_matrix,
            tpm_matrix=self.tpm_matrix,
            cpm_matrix=self.cpm_matrix
        )


class CalculateTPM:
    def __init__(
        self,
        sample_metrics: dict[str, SampleMetrics],
        filter_settings: FilterSettings,
        prep_method: PreparationMethod,
        context_name: str,
        output_directory: Path,
        per_million_scaling: int,
    ):
        self._sample_metrics: dict[str, SampleMetrics] = sample_metrics
        self._filter_settings: FilterSettings = filter_settings
        self._prep_method: PreparationMethod = prep_method
        self._context_name: str = context_name
        self._output_directory: Path = output_directory
        self._per_million_scaling: int = per_million_scaling
        self.__functions: list[Callable] = []
        
        self.__post_init__()
    
    def __call__(self, x):
        for function in self.__functions:
            result = function(x)
            if result is None or not result:
                return False
            return True
    
    def __post_init__(self):
        """
        These functions will be called after the object is initialized
        """
        self._output_directory.mkdir(parents=True, exist_ok=True)
        
        self._calculate_tpm()
        self._calculate_high_confidence()
    
    def _calculate_tpm(self):
        """
        This function will calculate the TPM for each sample
        """
        for key in self._sample_metrics.keys():
            current_metric: SampleMetrics = self._sample_metrics[key]
            study_number: str = current_metric.study_number
            entrez_ids: pd.Series = current_metric.entrez_ids
            size: pd.Series = current_metric.gene_size
            count_matrix: pd.DataFrame = current_metric.count_matrix
            gene_size: pd.Series = current_metric.gene_size
            
            tpm_matrix: pd.DataFrame = count_matrix.divide(gene_size, axis=0) * self._per_million_scaling
            current_metric.tpm_matrix = tpm_matrix
            current_metric.tpm_matrix.columns = count_matrix.columns
            
            # Write TPM matrix to file
            output_file: Path = self._output_directory / f"TPM_Matrix_{self._context_name}_{self._prep_method.value}_{study_number}_tpm.csv"
            output_df: pd.DataFrame = pd.concat([entrez_ids, current_metric.tpm_matrix], axis=1)
            output_df.to_csv(output_file)
            
            self._sample_metrics[key] = current_metric
    
    def _calculate_high_confidence(self):
        """
        This function will calculate the high confidence genes of each sample
        """
        N_exp = self._filter_settings.replicate_ratio
        N_top = self._filter_settings.replicate_ratio_high
        quantile = self._filter_settings.quantile
        
        for key in self._sample_metrics.keys():
            current_metric: SampleMetrics = self._sample_metrics[key]
            
            min_samples = np.round(N_exp * len(current_metric.tpm_matrix.columns))
            top_samples = np.round(N_top * len(current_metric.tpm_matrix.columns))
            test_boolean = pd.DataFrame(data={"gene": current_metric.entrez_ids})
            
            for i in range(len(current_metric.tpm_matrix.columns)):
                tpm_quantile = current_metric.tpm_matrix[:, i]
                tpm_quantile = tpm_quantile[tpm_quantile > 0]
                quantile_cutoff = np.quantile(tpm_quantile, q=1 - quantile / 100)
                test_boolean = pd.concat([
                    test_boolean,
                    pd.Series(np.where(current_metric.tpm_matrix.iloc[:, i] > quantile_cutoff, 1, 0))
                ], axis=1)
            
            test_boolean["gene"] = None
            filter_func = self._kOverA(k=min_samples, A=0.9)
            function_list = self._filter_fun(filter_func)
            keep_genes = self._gene_filter(test_boolean, [function_list])
            
            current_metric.entrez_ids = current_metric.entrez_ids[keep_genes]
            current_metric.gene_size = current_metric.gene_size[keep_genes]
            current_metric.count_matrix = current_metric.count_matrix[keep_genes, :]
            current_metric.tpm_matrix = current_metric.tpm_matrix[keep_genes, :]
            
            top_filter_func = self._kOverA(k=top_samples, A=0.9)
            top_function_list = self._filter_fun(top_filter_func)
            keep_top = self._gene_filter(test_boolean, [top_function_list])
            current_metric.entrez_ids_high_confidence = current_metric.entrez_ids[keep_top]
            
            self._sample_metrics[key] = current_metric
    
    def _kOverA(self, k, A, na_rm=True):
        """
        This is a reimplementation of R's genefilter::kOverA function
        
        R code: https://github.com/Bioconductor/genefilter/blob/e0b7d1614883990f3cc6f7fca071a814b1899d18/R/all.R#L5-L11
        ```
        """
        
        def filter_function(x):
            if na_rm:
                x = x[~np.isnan(x)]
            return np.sum(x > A) >= k
        
        return filter_function
    
    def _filter_fun(self, *args):
        """
        This is a reimplementation of R's genefilter::filterfun function. It uses a helper class, TPMFilterFun (above)
        
        R code: https://github.com/Bioconductor/genefilter/blob/e0b7d1614883990f3cc6f7fca071a814b1899d18/R/all.R#L133-L148
        """
        if len(args) == 1 and isinstance(args[0], list):
            function_list = args[0]
        else:
            function_list = list(args)
        self.__functions = function_list
        return self(*function_list)
    
    def _gene_filter(
        self,
        expression: Union[pd.DataFrame, np.ndarray],
        function_list: list[Callable]
    ):
        """
        This is a reimplementation of R's genefilter::genefilter function
        
        R code: https://github.com/Bioconductor/genefilter/blob/e0b7d1614883990f3cc6f7fca071a814b1899d18/R/all.R#L126-L131
        """
        if isinstance(expression, pd.DataFrame):
            expression = expression.values
        
        return np.apply_along_axis(function_list, axis=1, arr=expression)
    
    @property
    def sample_metrics(self):
        return self._sample_metrics


def create_sample_metrics(
    counts_matrix_file: Path,
    config_file: Path,
    info_file: Path,
    context_name: str
) -> dict[str, SampleMetrics]:
    config_object: pd.DataFrame = pd.read_excel(config_file, sheet_name=context_name)
    counts_matrix: pd.DataFrame = pd.read_csv(counts_matrix_file)
    counts_matrix = counts_matrix.sort_values("genes")
    
    gene_info: pd.DataFrame = pd.read_csv(info_file)
    gene_info["size"] = gene_info["end_position"] - gene_info["start_position"]  # add 'size' column
    gene_info = gene_info.sort_values("ensembl_gene_id")  # arrange by ensembl_gene_id
    gene_info = gene_info[gene_info["ensembl_gene_id"].isin(counts_matrix["genes"])]  # filter based on matching genes
    
    counts_matrix = counts_matrix.loc[gene_info["entrezgene_id"] != "-"]  # remove unnamed genes
    gene_info = gene_info[gene_info["entrezgene_id"] != "-"]  # remove unnamed genes
    genes = gene_info["entrezgene_id"]  # get gene names
    
    # Remove version numbers from ensembl id
    for i in range(len(genes)):
        row = genes.iloc[i]
        if re.search(r"\.", row):
            modified = row.split(".")[0]
            genes[i] = modified
    
    sample_metrics: dict[str, SampleMetrics] = {}
    groups = config_object["Group"].unique()
    for group in groups:
        metrics = SampleMetrics(count_matrix=genes)
        sample_metrics[group] = metrics
    
    # Add to group count matrices and insert lists
    for i in range(len(config_object["SampleName"])):
        entry = config_object["SampleName"][i]  # CELL-TYPE_S##R##
        group = config_object["Group"][i]
        
        if entry in counts_matrix.columns:
            # These values will be added to the sample_matrix object
            new_sample_matrix = counts_matrix[entry]
            new_fragment_length = config_object["FragmentLength"][i]
            new_layout = config_object["Layout"][i]
            new_sample_name = entry
            
            # Incase insert_size is None, set it to 0
            if new_fragment_length is None:
                new_fragment_length = 0
            
            # Update values
            sample_metrics[group].layout.append(new_layout)
            sample_metrics[group].sample_names.append(new_sample_name)
            sample_metrics[group].fragment_lengths.append(new_fragment_length)
            sample_metrics[group].count_matrix = pd.concat(
                [
                    sample_metrics[group].count_matrix,
                    new_sample_matrix
                ],
                axis=1
            )
        else:
            print(f"WARNING: {entry} not found in count matrix {str(counts_matrix_file)}")
    
    for group in groups:
        sample_matrix = sample_metrics[group].count_matrix
        sample_matrix = sample_matrix.iloc[:, 1:]
        sample_matrix = sample_matrix.apply(pd.to_numeric)
        sample_matrix.columns = sample_metrics[group].sample_names  # set column names to sample names
        
        sample_metrics[group].count_matrix = sample_matrix  # update counts matrix
        sample_metrics[group].num_samples = sample_matrix.shape[1]  # set number of samples
        sample_metrics[group].entrez_ids = gene_info["entrezgene_id"].astype(str)  # store entrez ids
        sample_metrics[group].gene_size = gene_info["size"]  # store gene size
        sample_metrics[group].study_number = group
    
    return sample_metrics


def calculate_paired_fpkm(
    count_column: np.ndarray,
    gene_size: np.ndarray,
    fragment_length: int
):
    effective_length = gene_size - fragment_length + 1
    N = count_column.sum()
    return np.exp(
        np.log(count_column)
        + np.lg(1e9)
        - np.log(effective_length)
        - np.log(N)
    )


def calculate_single_fpkm(
    count_column: np.ndarray,
    gene_size: np.ndarray
):
    count_column += 1
    gene_size += 1
    rate = np.log(count_column) - np.log(gene_size)
    return np.exp(
        rate
        - np.log(count_column.sum())
        + np.log(1e9)
    )


def calculate_fpkm(sample_metrics: dict[str, SampleMetrics]):
    for i, (key, metrics) in enumerate(sample_metrics.items()):
        layout = metrics.layout[i]
        assert layout in ["paired-end", "single-end"], \
            f"Invalid layout specified: {layout}. " \
            f"Must be one of ['paired-end', 'single-end']"
        
        count_matrix = metrics.count_matrix
        gene_size = metrics.gene_size
        if layout == "paired-end":
            fragment_lengths = metrics.fragment_lengths
            fpkm_matrix = np.column_stack(
                [
                    calculate_paired_fpkm(count_matrix[column].values, gene_size.values, fragment_lengths[i])
                    for column in count_matrix.columns
                ]
            )
            fpkm_matrix[np.isnan(fpkm_matrix)] = 0
            fpkm_matrix = pd.DataFrame(fpkm_matrix, columns=count_matrix.columns)
            metrics.fpkm_matrix = fpkm_matrix
        
        elif layout == "single-end":
            rpkm_matrix = np.column_stack(
                [
                    # Get column values
                    calculate_single_fpkm(count_matrix[column].values, gene_size.values)
                    for column in count_matrix.columns
                ]
            )
            rpkm_matrix[np.isnan(rpkm_matrix)] = 0
            rpkm_matrix = pd.DataFrame(rpkm_matrix, columns=count_matrix.columns)
            metrics.fpkm_matrix = rpkm_matrix
        print(metrics.fpkm_matrix)
    return sample_metrics


def calculate_z_score(
    sample_metrics: dict[str, SampleMetrics],
    normalization_technique: NormalizationTechnique,
) -> dict[str, SampleMetrics]:
    for key in sample_metrics.keys():
        current_metric: SampleMetrics = sample_metrics[key]
        if normalization_technique.value == "cpm":
            matrix = current_metric.cpm_matrix
        elif normalization_technique.value == "tpm":
            matrix = current_metric.tpm_matrix
        
        # Create a matrix with the same dimensions as `matrix`
        z_matrix = np.zeros(matrix.shape)
        
        for i in range(len(matrix.columns)):
            t_vector = matrix.iloc[:, i]
            log_vector = np.log2(t_vector)
            log_vector[np.isinf(log_vector)] = np.nan
            z_vector = sklearn.preprocessing.scale(log_vector, with_mean=True, with_std=True)
            z_matrix[:, i] = z_vector
        
        z_matrix = pd.DataFrame(z_matrix, columns=matrix.columns)
        
        current_metric.z_matrix = z_matrix
        sample_metrics[key] = current_metric
    return sample_metrics


def cpm_filter(
    sample_metrics: dict[str, SampleMetrics],
    filter_settings: FilterSettings,
    context_name: str,
    prep_method: PreparationMethod,
) -> dict[str, SampleMetrics]:
    """
    This is a python implementation R's edgeR::cpm function
    It will calculate RPKM values (as opposed to true CPM) for each gene in each sample
    This is done because TPM is more robust than vanilla CPM
    """
    N_exp = filter_settings.replicate_ratio
    N_top = filter_settings.replicate_ratio_high
    min_count = filter_settings.min_count
    
    for key in sample_metrics.keys():
        current_metric = sample_metrics[key]
        study_number = current_metric.study_number
        counts = current_metric.count_matrix
        entrez_ids = current_metric.entrez_ids
        gene_size = current_metric.gene_size
        library_size = counts.sum(axis=0)  # sum across rows


def tpm_filter(
    sample_metrics: dict[str, SampleMetrics],
    filter_settings: FilterSettings,
    context_name: str,
    prep_method: PreparationMethod,
    output_directory: Path,
    per_million_scaling: int = 1e6
) -> dict[str, SampleMetrics]:
    return CalculateTPM(
        sample_metrics=sample_metrics,
        filter_settings=filter_settings,
        prep_method=prep_method,
        context_name=context_name,
        output_directory=output_directory,
        per_million_scaling=per_million_scaling
    ).sample_metrics


def zfpkm_filter():
    pass


def umi_filter():
    pass


def filter_counts(
    sample_metrics: dict[str, SampleMetrics],
    filter_settings: FilterSettings,
    context_name: str,
    prep_method: PreparationMethod,
    output_directory: Path
) -> dict[str, SampleMetrics]:
    if filter_settings.technique == "cpm":
        return cpm_filter()
    elif filter_settings.technique == "zfpkm":
        return zfpkm_filter()
    elif filter_settings.technique == "quantile":
        return tpm_filter(
            sample_metrics=sample_metrics,
            filter_settings=filter_settings,
            context_name=context_name,
            prep_method=prep_method,
            output_directory=output_directory
        )
    elif filter_settings.technique == "umi":
        return umi_filter()
    else:
        raise ValueError(
            f"Invalid technique: {filter_settings.technique}. Must be one of ['cpm', 'zfpkm', 'quantile', 'umi']")


def main(
    counts_matrix_file: Path,
    config_file: Path,
    out_file: Path,
    info_file: Path,
    context_name: str,
    prep: str = "total",
    replicate_ratio: float = 0.5,
    batch_ratio: float = 0.5,
    replicate_ratio_high: float = 0.9,
    batch_ratio_high: float = 0.9,
    technique: str = "quantile",
    quantile: float = 0.9,
    min_count: int = 10,
    min_zfpkm: int = -3
):
    # Condense filtering options
    filter_settings = FilterSettings(
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        quantile=quantile,
        min_count=min_count,
        min_zfpkm=min_zfpkm,
        replicate_ratio_high=replicate_ratio_high,
        batch_ratio_high=batch_ratio_high,
        technique=technique,
    )
    
    prep: PreparationMethod = PreparationMethod(prep)
    if prep.value.lower() == "scrna":
        filter_settings.technique = "umi"
        print("Note: Single cell filtration does not normalize and assumes counts are counted with UMI")
    
    sample_metrics = create_sample_metrics(
        counts_matrix_file=counts_matrix_file,
        config_file=config_file,
        info_file=info_file,
        context_name=context_name,
    )
    # calculate_fpkm(sample_metrics)
    
    result = filter_counts(
        sample_metrics=sample_metrics,
        filter_settings=filter_settings,
        context_name=context_name,
        prep_method=prep,
    )
    print(result["S3"].count_matrix)
    
    # expressed_genes: list = []
    # top_genes: list = []
    # for key in sample_metrics.keys():
    #     expressed_genes.append(sample_metrics[i]["Entrez"])
    #     top_genes.append(sample_metrics[i]["Entrez_hc"])


if __name__ == "__main__":
    main(
        counts_matrix_file=Path("/Users/joshl/Downloads/COMO_python/input/gene_counts_matrix_total_naiveB.csv"),
        info_file=Path("/Users/joshl/Downloads/COMO_python/input/gene_info.csv"),
        config_file=Path("/Users/joshl/Downloads/COMO_python/input/trnaseq_data_inputs_auto.xlsx"),
        out_file=Path("/Users/joshl/Downloads/COMO_python/output/naiveB.csv"),
        context_name="naiveB",
        technique="quantile"
    )
