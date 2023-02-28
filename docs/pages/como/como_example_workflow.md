---
title: Example Workflow
sidebar: como_sidebar
permalink: como_example_workflow.html
folder: como
last_updated: September 14, 2022
---

## Overview
This example will identify drug targets for [Rheumatoid Arthritis](https://en.wikipedia.org/wiki/Rheumatoid_arthritis) using GSMN of naïve [CD4+ T-cell subtypes](https://en.wikipedia.org/wiki/T_helper_cell) 

Follow the steps found at [Starting the Container](como_getting_started.html#starting-the-container) to start the container.

## Step 1

The outputs of [STAR aligner](https://github.com/alexdobin/STAR) with the `--GeneCounts` argument can be interfaced directly with COMO. Under the `{{ site.data.terms.data_dir }}/STAR_output/` folder, we provide a template structure with naïve B control cells from bulk RNA-sequencing experiments found in [NCBI's Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). We can run the file `{{ site.data.terms.py_dir }}/rnaseq_preprocess.py` with the `create_count_matrix` argument set to `True` to merge the counts from each replicate in each study into a single matrix.<br>This follows **Step 1** in the COMO container.

Alternatively, you can provide a count matrix directly and run with `create_count_matrix` set to `False` to just fetch the necessary information required for normalization

## Step 2

In the `{{ site.data.terms.config_sheets }}` folder, we include the following files:
- `microarracy_data_inputs.xlsx` contains GEO accession numbers of microarray data
- `proteomics_data_inputs.xlsx` contains sample names of proteomics data

Protein abundance is located under `{{ site.data.terms.config_sheets }}/ProteomicsDataMatrix_Naive.csv`. Sample names for bulk RNA-sequencing data is given in `{{ site.data.terms.config_sheets }}/bulk_data_inputs.xlsx`. Running `{{ site.data.terms.py_dir }}/rnaseq_preprocess.py` will create `{{ site.data.terms.results_dir }}/Gene_Info_naiveB.csv` and `{{ site.data.terms.data_matrices_dir }}/BulkRNAseqMatrix_naiveB.csv` if `create_counts_matrix` is set to `True`.

`{{ site.data.terms.config_sheets }}/microarray_data_inputs.xlsx` includes microarray samples of naice CD4+ T cells from 
[GSE22886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22886), 
[GSE43005](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43005), 
[GSE22045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22045), and 
[GSE24634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24634).
The file `{{ site.data.terms.config_dir }}/proteomics_data_inputs.xlsx` contains sample names of naive CD4+ T cells, with its results found in the file `{{ site.data.terms.data_matrices_dir }}/naiveB/protein_abundance_naiveB.csv`

Using `{{ site.data.terms.py_dir }}/merge_xomics.py`, you can specify any number of available data sources (microarray, bulk RNA-seq, and proteomics) as inputs. You can also set the `expression_requirement` parameter, which defines the minimum number of data sources that must have gene expression above the threshold limit for the said gene to be considered active

{% include note.html content="If a gene is not supported by a data source or platform, the `expression_requirement` value will decrease by one for each input data source that does not support the gene." %}

Running **Step 1** in the COMO container will generate "gene activity" files based on transcriptomics and proteomics data, as described by <a href="https://doi.org/10.1038/s41540-020-00165-3" target="_blank"><cite>Puniya et al., 2020</cite></a>

This will save final output in `GeneExpression_Naive_merged.csv` and its path in `step1_results_files.json`

## Step 3
Our pipeline includes a modified version of the [Recon3D](https://doi.org/10.1038/nbt.4072) model to use as a reference for model contextualization. The modified version of Recon3D is available at `{{ site.data.terms.data_dir }}/GeneralModelUpdatedV2.mat`.

**Step 4** in COMO will use `GeneExpression_Naive_merged.csv` (from **Step 1**, above) in combination with `{{ site.data.terms.data_dir }}/GeneralModelUpdatedV2.mat` to construct a cell-type specific model of Naive CD4+ cells.<br>
We can use this model in the next steps. However, we advise users to properly investigate, manually curate, and reupload the refined version in `{{ site.data.terms.results_dir }}` to use in **Step 4**.

## Step 4
We used a dataset ([GSE56649](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56649)) of Rheumatoid Arthritis to identify differentially expressed genes (disease genes). We defined accession IDs of this dataset in the inputs found under `{{ site.data.config_dir }}/naiveB/disease/gene_counts_matrix_arthritis_naiveB.csv`.

This step will generate files `{{ site.data.terms.results_dir }}/naiveB/Disease_UP_GSE56649.txt` and `{{ site.data.terms.results_dir }}/naiveB/Disease_DOWN_GSE56649.txt`, and their paths, at `step2_results_files.json`. Finally, this step will create a `disease_files` variable that will include aths of files for up- and down- regulated genes.

## Step 5
This step will use the model (constructed in **Step 2**/uploaded curated version) and perform knock-out simulations of genes overlapping with the drug-target data file obtained from the [ConnectivityMap](https://www.broadinstitute.org/connectivity-map-cmap) database. We refined the drug target-data file and included it at `{{ site.data.terms.data_dir }}/RepurposingHub.txt`.

This step will use the following files:
- `Naive_SpecificModel.json` (or a pre-curated version uploaded as `NaiveModel.mat`)
- `Disease_UP_GSE56649.txt` 
- `Disease_DOWN_GSE56649.txt` 
- `RepurposingHub.txt`

The final output files will include drug targets ranked based on Perburbation Effect Score (PES) as described by <a href="https://doi.org/10.1038/s41540-020-00165-3" target="_blank"><cite>Puniya et al., 2020</cite></a>.

The output file `{{ site.data.terms.data_dir}}/output/d_score.csv` will contain Entrez IDs of ranked genes and their corresponding PES. The file `drug_score.csv` will contain PES ranked drug targets (Entrez IDs and Gene Symbols) with mapped repurposed drugs.
