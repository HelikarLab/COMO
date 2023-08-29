---
title: Getting Started
sidebar: como_sidebar
permalink: index.html
summary: This is an overview of how to get started with the COMO container
last_updated: Sept 15, 2022
---

## What is COMO?

COMO stands for **C**onstraint-based **O**ptomization of **M**etabolic **O**bjectives

This is a [Jupyter Notebook](https://jupyter.org/)-based pipeline to build
context-specific, [constraint-based metabolic models (CBMM)](https://en.wikipedia.org/wiki/Metabolic_network_modelling)
from a single source, or a combination of sources, from the following "-omics" data:

- [Bulk RNA-sequencing](https://en.wikipedia.org/wiki/RNA-Seq)
- [Single-Cell](https://en.wikipedia.org/wiki/Single_cell_sequencing) RNA-sequencing
- [Mass Spectrometry](https://en.wikipedia.org/wiki/Mass_spectrometry) [Proteomics](https://en.wikipedia.org/wiki/Proteomics)
- [Microarray](https://en.wikipedia.org/wiki/Microarray)

Ultimately, COMO serves as a platform to build constraint-based models, identify drug targets, and discover potentially
repurposable drugs for metabolism-impacting diseases.

COMO does not require any amount of programming experience to create models. However, every step of the pipeline is
packaged in its own python file to promote accessible modification, addition, or replacement of analysis steps. The
Jupyterlab container comes pre-loaded with the most popular [R](https://www.r-project.org/)
and [Python](https://www.python.org/) libraries, but if you would like to use a library and cannot install it,
please [open a new issue](https://github.com/HelikarLab/COMO/issues)
on [our Github page](https://github.com/HelikarLab)!

In this pipeline, the term "context" refers to a specific state that can be subset from a genome-scale metabolic model
using experimental data. This context can be a cell or tissue type in a specific experimental state, such as a
"healthy" or "diseased" state.

For drug perturbation scoring of a specific cell or tissue type, it is only necessary to build a constraint-based
metabolic model (Steps 1 and 2) for the healthy control state. Differential gene expression (Step 3) will be used to for
disease analysis so that multiple diseaes can be analyzed using the same constraint-based model.

## Before You Start

The proper data must be provided, dependent on what analysis you would like to run

### RNA-Sequencing

- A folder named "COMO_input" in the `{{ site.data.terms.data_dir }}` directory is available. Proper inputs can be
  generated using our [SnakeMake Pipeline](https://github.com/HelikarLab/FastqToGeneCounts), designed specifically for
  COMO.
- If processing RNA-sequencing data with an alternate procedure, or importing a pre-made gene count matrix, follow the
  instructions found in **Step 1** of the COMO container.

{% include note.html content="RNA-sequencing data can be [bulk](https://en.wikipedia.org/wiki/RNA-Seq)
or [single-cell](https://en.wikipedia.org/wiki/Single_cell_sequencing) sequencing, but our SnakeMake pipeline is
currently only available for bulk RNA-sequencing." %}

### Proteomics

- A matrix of measurements where the rows are proteins in Entrez format, and columns are sample names

{% include note.html content="A proteomics workflow is included in the COMO container." %}

### Microarray

{% include warning.html content="Microarray has become mostly obsolete. Another data source should be used if possible,
but is included in COMO if you'd like to use it." %}

- Results must be uploaded to [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
- The only thing required by COMO is a configuration file with: GSE, GSM, and GPL codes
- A template can be found
  at: [`{{ site.data.terms.data_dir }}/config_sheets/microarray_data_inputs_template.xlsx`](https://github.com/HelikarLab/COMO/blob/master/main/data/config_sheets/microarray_data_inputs_template.xlsx)

## Configuration Information

Configuration Excel files should be uploaded
to [`{{ site.data.terms.data_dir }}/config_sheets/`](https://github.com/HelikarLab/COMO/tree/master/main/data/config_sheets).
The sheet names in these configuration files should correspond to the context (tissue name, cell name, control, etc.)
where each sheet contains the sample names to include in that context-specific model. These sample names should
correspond to the sample (column names in the source data matrix, which should be uploaded (or output)
in [`{{ site.data.terms.data_dir }}/data_matrices/<model name>`](https://github.com/HelikarLab/COMO/tree/master/main/data/data_matrices)
of the COMO container.

{% include note.html content="As stated in [this issue](https://github.com/HelikarLab/COMO/issues/46), we will be moving
to use pure CSV files for ease of use in COMO. There is no timeline for this, as of now."%}

In the Docker image, some exemplary input files are included to build metabolic models of naïve, Th1, Th2, and Th17
subtypes, and identify drug targets for Rheumatoid arthritis. You can follow the documentation in the container, the
format of these files, and the template files to create your own input files.

## Obtaining a Gurobi License

From
the [MathWorks Website](https://www.mathworks.com/products/connections/product_detail/gurobi-optimizer.html#:~:text=The%20Gurobi%20Optimizer%20is%20a,implementations%20of%20the%20latest%20algorithms.):
> The Gurobi Optimizer is a state-of-the-art solver for mathematical programming. THe solvers in the Gurobi Optimizer
> were designed from the ground up to exploit modern architectures and multicore processors, using the most advanced
> implementations of the latest algorithms

This is to say, Gurobi will help when attempting to calculate linear programming solutions done during flux balance
analysis calculations

A Gurobi license is free for academic use. To obtain a license, perform the following:

1. Access the [Web License Manager](https://license.gurobi.com/manager/keys), and log in (or create an account)
2. Click on "API Keys" on the left sidebar
3. Click "Create API Key"
4. Enter the relevant information. The "Application Name" can be anything; I have mine set to "COMO", with an empty
   description
5. Click "Create"
6. You **MUST** download the license file now, as it is not available later. Creating new licenses is easily done,
   however.
7. Save this file to a location on your computer, and note the location. This will be used in the next step.

## Accessing COMO

Navigate to [http://localhost:8888](http://localhost:8888) in your browser

## Working in the Notebook

Open the Jupyter Notebook file, found at `COMO.ipynb` in the web interface.

Configuration files should be uploaded to `{{ site.data.terms.data_dir }}/config_files` and data files
to `{{ site.data.terms.data_dir }}/data_matrices`, according to instructions in the notebook and provided templates

Update the file names in the Jupyter Notebook accordingly
