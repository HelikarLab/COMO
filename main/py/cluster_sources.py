from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# read and translate R functions
r_file_source = open("rscripts/cluster_sources.R", "r").read()

cluster_sources = SignatureTranslatedAnonymousPackage(r_file_source, "cluster_sources")

results_dir = "/Users/joshl/docker/madrid/local_files/results"
context_names = ["immNK", "naiveB"]
source_type = "zFPKM"
use_trna = True
use_mrna = True
binarize_data = True

cluster_sources.cluster_sources_main(
    results_dir=results_dir,
    context_names=context_names,
    source_type=source_type,
    use_trna=use_trna,
    use_mrna=use_mrna,
    binarize_data=binarize_data
)
