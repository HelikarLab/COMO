"""
This file successfully creates cluster names using scanpy's 3kPBMC tutorial
"""

"""
This file is the beginning of imlementing a gradient boosting classification algorithm using scikit-learn.
"""
from pathlib import Path
import scanpy as sc
from scanpy import AnnData


def get_annotated_filtered_dataframe(data_dir: Path | str, debug: bool = False) -> sc.AnnData:
    if isinstance(data_dir, str):
        data_dir = Path(data_dir)
    
    # ----- Reading -----
    if debug:
        print("Reading")
    adata: AnnData = sc.read_10x_mtx(data_dir, var_names="gene_symbols", cache=True)
    
    # ----- Filtering -----
    if debug:
        print("Filtering")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Annotate the group of mitochondrial genes as "mt"
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    
    # Remove genes with counts over 2500
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    # Remove genes with percentage of mitochondrial counts over 5%
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Calculate log1p
    sc.pp.log1p(adata)
    
    # Calculate highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Set pre-filtered data as raw
    adata.raw = adata
    
    # Set highly variable genes as filtered data
    adata = adata[:, adata.var.highly_variable]
    
    # ----- Normalization -----
    if debug:
        print("Normalizing")
    # Remove sources of variation
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)
    
    # ----- Non-Annotated Clustering -----
    if debug:
        print("Clustering")
    sc.tl.pca(adata, svd_solver="arpack")
    
    # Create clusters
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata)
    cluster_names: list[str] = [str(i) for i in range(adata.obs["leiden"].nunique())]
    adata.rename_categories("leiden", cluster_names)
    
    # ----- Plotting -----
    if debug:
        print("Plotting")
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata, init_pos="paga")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    
    return adata


def main():
    # Set scanpy verbosity to minimum
    sc.settings.verbosity = 0
    data_dir: Path = Path("/Users/joshl/PycharmProjects/MADRID/main/py/single_cell/single_cell_data/3kPBMC")
    adata: AnnData = get_annotated_filtered_dataframe(data_dir, debug=True)
    
    # plot = sc.pl.umap(adata, color=["leiden", "CD19"], use_raw=False, show=False)
    
    adata.rename_categories("leiden", ['CD4 T', 'CD14 Monocytes', 'B', 'CD8 T', 'NK', 'FCGR3A Monocytes', 'Dendritic',
                                       'Megakaryocytes'])
    sc.pl.umap(adata, color="leiden", title="", legend_loc="on data", frameon=False, show=True)
    
    print("DONE")


if __name__ == '__main__':
    main()
