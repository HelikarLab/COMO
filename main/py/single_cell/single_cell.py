import pandas as pd
import scanpy as sc
import gzip
from scanpy import pp as preprocess
from scanpy import pl as plot
from scanpy import tl as tools
from scanpy import AnnData
from pathlib import Path
from typing import Literal, Iterable
import pickle

b_cell_markers: dict[str, int] = {
    "CD19": -1, "CD20": -1, "CD21": -1, "CD22": -1,
    "CD23": -1, "CD24": -1, "CD27": -1, "CD79a": -1
}
t_cell_markers: dict[str, int] = {
    "CD3": -1, "CD4": -1, "CD8": -1, "CD45": -1,
    "CD69": -1, "CD25": -1, "CD127": -1, "CD279": -1
}


def filter_data(adata: AnnData) -> AnnData:
    # Perform filtering on data
    preprocess.filter_cells(adata, min_genes=200, inplace=True)
    preprocess.filter_genes(adata, min_cells=0, inplace=True)
    return adata


def preprocess_data(adata: AnnData, input_path: Path | str, force_preprocess: bool = False) -> AnnData:
    # Perform preprocessing
    # Normalized total
    preprocessed_data: Path = Path("./cache", f"data-{input_path.stem}-preprocess.pickle")
    
    if not preprocessed_data.is_file() or force_preprocess:

        tools.pca(adata, svd_solver='arpack')
        preprocess.neighbors(adata, n_neighbors=10, n_pcs=40)
        tools.leiden(adata, key_added="leiden", neighbors_key="neighbors")
        
        tools.paga(adata)
        sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
        sc.tl.umap(adata, init_pos='paga')

        tools.rank_genes_groups(adata, "leiden", method="wilcoxon")
        
        plot.rank_genes_groups(adata, n_genes=25, sharey=False)

        # preprocess.normalize_total(adata, inplace=True)
        # preprocess.log1p(adata)
        # preprocess.highly_variable_genes(adata)
        # preprocess.neighbors(adata, use_rep="X", )
        # tools.leiden(adata)
        # tools.rank_genes_groups(adata, "leiden")
        # tools.pca(adata)
        
        bytes_ = pickle.dumps(adata)
        with open(preprocessed_data, "wb") as o_stream:
            o_stream.write(bytes_)
    else:
        print("\tLoading normalized data from cache")
        # adata = pickle.loads(normalized_total_path.read_bytes())
        with open(preprocessed_data, "rb") as i_stream:
            adata = pickle.load(i_stream)
    
    return adata


def plot_umap(adata: AnnData, color_: str = "leiden", legend_loc: str = "right margin", show_plot: bool = True):
    tools.umap(adata, init_pos="paga")
    umap = plot.umap(adata, color=color_, legend_loc=legend_loc, show=show_plot, title="UMAP")
    return umap


def plot_tsne(adata: AnnData, color_: str = "leiden", legend_loc: str = "right margin", show_plot: bool = True):
    tools.tsne(adata, n_pcs=2)
    tsne = plot.tsne(adata, color=color_, show=show_plot, legend_loc=legend_loc, title="tSNE")
    return tsne


def cluster(adata: AnnData, method: Literal["umap", "tsne", "both"], show_plot: bool = False) -> AnnData:
    # Cluster using the specified method
    tools.paga(adata)
    plot.paga(adata, plot=False)  # Only compute the layout
    cluster_names: list[str] = [str(i) for i in range(adata.obs["leiden"].nunique())]
    adata.rename_categories("leiden", cluster_names)
    
    match method:
        case "umap":
            plot_umap(adata=adata, color_="leiden", show_plot=show_plot)
        case "tsne":
            plot_tsne(adata=adata, color_="leiden", show_plot=show_plot)
        case "both":
            print("\tPlotting UMAP")
            plot_umap(adata=adata, color_="leiden", show_plot=show_plot)
            print("\tPlotting tSNE")
            plot_tsne(adata=adata, color_="leiden", show_plot=show_plot)
    
    return adata


def get_ranked_genes(adata: AnnData, group: str | Iterable[str] = "0") -> pd.DataFrame:
    results: pd.DataFrame = sc.get.rank_genes_groups_df(adata, group=group)
    return results


def read_data(path: Path | str) -> AnnData:
    if isinstance(path, str):
        read_directory: Path = Path(path)
    else:
        read_directory: Path = path
    
    if not read_directory.is_dir():
        raise ValueError(f"Variable `path` must be a directory")
    elif not read_directory.exists():
        raise ValueError(f"Unable to find a directory named {str(read_directory)}")
    
    try:
        adata: AnnData = sc.read_10x_mtx(
            read_directory,  # the directory with the `.mtx` file
            var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
            cache=True)  # write a cache file for faster subsequent reading
    
        adata.var_names_make_unique()
    except FileNotFoundError:
        print(f"gzipped matrix file not found, gzipping the matrix.mtx file under {str(read_directory)}")
        
        # Iterate through each file in the read_directory:
        for file in read_directory.glob("*"):
            if file.name in ["barcodes.tsv", "features.tsv", "matrix.mtx"]:
                gzip_write_filepath = Path(read_directory, f"{file.name}.gz")
                with open(file, "rb") as reader, gzip.open(gzip_write_filepath, "wb") as writer:
                    writer.writelines(reader)
        
        adata: AnnData = read_data(read_directory)
    
    return adata


def pbmc_example():
    path_: Path = Path("data/GSE128243 NK Cell/Unstimulated 3")
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    adata: AnnData = read_data(path_)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    # Testing below parameters
    tools.pca(adata)
    preprocess.neighbors(adata, n_neighbors=10, n_pcs=40)
    tools.leiden(adata, key_added="leiden")
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
    sc.tl.umap(adata, init_pos='paga')
    
    sc.pl.umap(adata, color=['leiden'], legend_loc="on data")

    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    adata.write(path_)


def main():
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(dpi=1000, vector_friendly=True, format="svg")
    
    # Path to folder
    input_path: Path = Path("data/GSE121267 T Cell Dataset/Donor 1")
    
    print("Reading")
    adata: AnnData = read_data(input_path)
    
    print("Filtering")
    # adata = filter_data(adata=adata)
    
    print("Preprocessing")
    adata = preprocess_data(adata=adata, input_path=input_path, force_preprocess=True)
    
    # print("Clustering and plotting")
    # adata = cluster(adata=adata, method="both", show_plot=True)
    
    # print("Collecting highly ranked genes")
    # groups: list[str] = [str(i) for i in adata.obs["leiden"].unique().tolist()]
    # ranked_genes: pd.DataFrame = sc.get.rank_genes_groups_df(adata=adata, group=groups)
    # ranked_genes.to_csv(input_path / "ranked_genes.csv")
    
    x = 0


if __name__ == '__main__':
    pbmc_example()
    # main()
