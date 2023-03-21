"""
This file is the beginning of imlementing a gradient boosting classification algorithm using scikit-learn.
"""

# import gzip
from pathlib import Path
from typing import Union

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

import scanpy as sc
from scanpy import AnnData


def get_annotated_filtered_dataframe(data: Union[Path, str], debug: bool = False) -> sc.AnnData:
    if isinstance(data, str):
        data = Path(data)

    # ----- Reading -----
    if debug:
        print("Reading")
    
    if data.is_file():
        if data.suffix == ".h5":
            adata: AnnData = sc.read_hdf(data, key="matrix")
    elif data.is_dir():
        # If a "KeyError" occurs here, unzip the files and read them again. The "features.tsv" file should be named "genes.tsv". Everything else stays the same
        adata: AnnData = sc.read_10x_mtx(data, cache=True)

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

    sc.tl.pca(adata)
#     print(adata.shape)
#     for i in range(1, 16):
#         print(f"Trying {i}")
#         try:
#             sc.tl.pca(adata, svd_solver="arpack", n_comps=i)
#         except ValueError:
#             pass

    # Create clusters
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=0.8, )
    cluster_names: list[str] = [str(i) for i in range(adata.obs["leiden"].nunique())]
    adata.rename_categories("leiden", cluster_names)

    # ----- Plotting -----
    if debug:
        print("Plotting")
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata, init_pos="paga")
    sc.tl.tsne(adata, learning_rate=200)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    
    return adata


def gradient_booster(adata: AnnData):
    import xgboost
    import sklearn
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import accuracy_score

    # Convert adata.obs["leiden"] to integers
    adata.obs["leiden"] = adata.obs["leiden"].astype(int)

    # ----- Gradient Boosting -----
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(adata.X, adata.obs["leiden"], test_size=0.2, random_state=42)

    # Create gradient boosting classifier
    classifier = xgboost.XGBClassifier(n_estimators=20, max_depth=3)

    # Train model
    classifier.fit(X_train, y_train)

    # Predict on test data
    y_pred = classifier.predict(X_test)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))


def main():
    print("Starting")
    # Set scanpy verbosity to minimum
    sc.settings.verbosity = 0
    
    # NCBI aligned data
    data_dir: Path = Path("/Users/joshl/PycharmProjects/MADRID/main/py/single_cell/single_cell_data/T Cells/GSE121267/GSM3430549")
    
    # self-aligned data
    # data_dir: Path = Path("/Users/joshl/PycharmProjects/MADRID/main/py/single_cell/single_cell_data/CD4 T Cells/GSE203217/GSM6164865")
    adata: AnnData = get_annotated_filtered_dataframe(data=data_dir, debug=False)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    
    sc.pl.umap(adata, color=["leiden", "CD4"], show=True)
    sc.pl.tsne(adata, color=["leiden", "CD4"], show=True)

    print("DONE")


if __name__ == '__main__':
    main()
