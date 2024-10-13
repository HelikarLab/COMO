from pathlib import Path

import rpy2_api

# read and translate R functions
r_file_path = Path("./rscripts/cluster_sources.R")
results_dir = "/Users/joshl/docker/madrid/local_files/results"
context_names = ["immNK", "naiveB"]
source_type = "zFPKM"
use_trna = True
use_mrna = True
binarize_data = True

cluster_sources = rpy2_api.Rpy2(
    r_file_path=r_file_path,
    results_dir=results_dir,
    context_names=context_names,
    source_type=source_type,
    use_trna=use_trna,
    use_mrna=use_mrna,
    binarize_data=binarize_data
).call_function("cluster_sources_main")
