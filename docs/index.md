---
title: Getting Started
sidebar: madrid_sidebar
tags: [getting_started]
permalink: index.html
summary: This is an overview of how to get started with the MADRID container
last_updated: Sept 15, 2022
---

## What is MADRID?
MADRID stands for "**M**et**A**bolic **D**rug **R**epurposing **ID**entification".

This is a [Jupyter Notebook](https://jupyter.org/)-based pipeline to build context-specific, [constraing-based metabolic models (CBMM)](https://en.wikipedia.org/wiki/Metabolic_network_modelling) from a single source, or a combination of sources, from the following "-omics" data:
- [Bulk RNA-sequencing](https://en.wikipedia.org/wiki/RNA-Seq)
- [Single-Cell](https://en.wikipedia.org/wiki/Single_cell_sequencing) RNA-sequencing
- [Mass Spectrometry](https://en.wikipedia.org/wiki/Mass_spectrometry) [Proteomics](https://en.wikipedia.org/wiki/Proteomics)
- [Microarray](https://en.wikipedia.org/wiki/Microarray)

Ultimately, this serves as a platform to use these models to identify drug targets and potentially repurposable drugs for metabolism-impacting diseases.

MADRID does not require any amount of programming experience to create models. However, every step of the pipeline is packaged in its own .py file to promote accessible modification, addition, or replacement of analysis steps. The Jupyterlab container comes pre-loaded with the most popular [R](https://www.r-project.org/) and [Python](https://www.python.org/) libraries, but if you would like to use a library and cannot install it, please [open a new issue](https://github.com/HelikarLab/MADRID/issues) on [our Github page](https://github.com/HelikarLab/MADRID)!

In this pipeline, the term "context" refers to a specific state that can be subset from a genome-scale metabolic model using experimental data. This context can be a specific cell or tissue type in a specific experimental state, such as a control or a disease.

For drug perturbation scoring of a specific cell-type or tissue, it is only necessary to build a CBMM (Steps 1 to 5) for the healthy control state. Differential gene expression (Step 6) will be used to for disease analysis so that multiple diseaes can be analyzed using the same CBMM.

## Before You Start

The proper data must be provided, dependent on what analysis you would like to run

### RNA-Sequencing
- A folder named "MADRID_input" in the `{{ site.data.terms.data_dir }}` directory is available. Proper inputs can be generated using our [SnakeMake Pipeline](https://github.com/HelikarLab/FastqToGeneCounts), designed specifically for MADRID.

{% include  note.html content="RNA-sequencing data can be [bulk](https://en.wikipedia.org/wiki/RNA-Seq) or [single-cell](https://en.wikipedia.org/wiki/Single_cell_sequencing) sequencing, but our SnakeMake pipeline is currently only available for bulk RNA-sequencing." %}

- If processing RNA-sequencing data with an alternate procedure, or importing a pre-made gene count matrix, follow the instructions found in **Step 1** of the MADRID container.

### Proteomics
- A matrix of measurement where rows are proteins in Entrez format, and columns are arbitrary sample names

{% include important.html content="A proteomics workflow is included in the MADRID container." %}

### Microarray
{% include warning.html content="Microarray has become mostly obsolete. Another data source should be used if possible, but is included in MADRID in case it is required." %}
- Results must be uploaded to [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
- The only thing required by MADRID is a configuration file with: GSE, GSM, and GPL codes
- A template can be found at: [`{{ site.data.terms.data_dir }}/config_sheets/microarray_data_inputs_template.xlsx`](https://github.com/HelikarLab/MADRID/blob/master/main/data/config_sheets/microarray_data_inputs_template.xlsx)

## Identifying Drug Targets
The following is a list of steps to identify drug targets. Stop after **Step 3** if you are building a context-specific model for other purposes

1. Preprocess Bulk RNA-sequencing data by converting gene counts from [STAR aligner](https://github.com/alexdobin/STAR) into a unified matrix, fetch necessary info about each gene required for normalization, and generate a configuration sheet.
2. Analyze any combination of Microarray, RNA-sequencing (total, polyA, or single-cell), or proteomics data, and output a list of active genes for each strategy/method.
3. Check for consensus amongst strategies according to desired rigor and merge into a singular set of active genes.
4. Create tissue specific models based on the list of active genes. If required, the user can manually refine these models and supply them in **Step 4** of the MADRID container.
5. Identify differential gene expressions from disease datasets using either microarray or bulk RNA-sequencing transcriptomics information.
6. Identify drug targets and repurposable drugs. This step consists of four substeps:
- Mapping drugs on automatically-created or user-supplied models
- Knock-out simulation
- Compare simulation results of perturbed and unperturbed models
- Integrate with disease genes and score drug targets

## Configuration Information
Configuration Excel files should be uploaded to [`{{ site.data.terms.data_dir }}/config_sheets/`](https://github.com/HelikarLab/MADRID/tree/master/main/data/config_sheets). The sheet names in these configuration files should correspond to the context (tissue name, cell name, control, etc.) where each sheet contains the sample names to include in that context-specific model. These sample names should correspond to the sample (column names in the source data matrix, which should be uploaded (or output) in [`{{ site.data.terms.data_dir }}/data_matrices/<model name>`](https://github.com/HelikarLab/MADRID/tree/master/main/data/data_matrices) of the MADRID container.

{% include note.html content="As stated in [this issue](https://github.com/HelikarLab/MADRID/issues/46), we will be moving to use pure CSV files for ease of use in MADRID. There is no timeline for this, as of now."%}

In the Docker image, some exemplary input files are included to build metabolic models of naïve, Th1, Th2, and Th17 subtypes, and identify drug targets for Rheumatoid arthritis. You can follow the documentation in the container, the format of these files, and the template files to create your own input files.

## Starting the Container

{% include warning.html content="If you terminate your session after running the Docker, any changes you make WILL NOT BE SAVED! Please mount a local directory to the Docker image as instructed on [Docker's website](https://docs.docker.com/storage/volumes/) to prevent data loss." %}

### Choosing a Tag
Several tags are available for use. [Tags](https://docs.docker.com/engine/reference/commandline/tag/) are Docker's method of versioning. For this package, several tags are available (as described below), but the most useful one is most likely `latest`.

A list of available tags can be found [here](https://github.com/HelikarLab/MADRID/pkgs/container/madrid/versions)

|      Tag      |                 Description                 |
|:-------------:|:-------------------------------------------:|
|   `latest`    |            Latest stable release            |
| `master-1.0`  |  Any stable changes under the 1.0 release   |

To start the container, [docker](https://docs.docker.com/) must be installed.

Execute the following tasks in a terminal to download the container, where `CHOSEN-TAG` is the tag chosen from the table above

```bash
sudo docker login
sudo docker pull ghcr.io/helikarlab/madrid:CHOSEN-TAG
```

For example, if I chose the `latest` tag, I would execute the following:
```bash
sudo docker login
sudo docker pull ghcr.io/helikarlab/madrid:latest
```

Now we can run the container.
```bash
docker run \
  --cpus=4 \
  -p 8888:8888 \
  --volume=$HOME/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
  --volume=$HOME/madrid_local:/home/jovyan/main/data/local_files \
  --name madrid \
  --rm \
  -it \
  "ghcr.io/helikarlab/madrid:CHOSEN-TAG"
```
{% include important.html content="Make sure to use the tag you chose from the table in the `docker run...` command, replacing `CHOSEN-TAG` with your tag." %}

Below is a list of options included in the above `docker run` command. These can be removed or changed as-needed, as long as you are aware of their implications
If you have questions, don't be afraid to ask on our [Issues Page](https://github.com/HelikarLab/MADRID/issues)!

|                 Option                 |  Required?  |                                                       Description                                                        |
|:--------------------------------------:|:-----------:|:------------------------------------------------------------------------------------------------------------------------:|
|               `--cpus=4`               |     No      |                   The number of CPUs to allocate to the container. This is optional, but recommended.                    |
|             `-p 8888:8888`             |     Yes     |                                        The port to use for the Jupyter Notebook.                                         |
|      `--volume=$HOME/gurobi.lic`       |     Yes     | The path to your Gurobi license file. This is required to run the [Gurobi Optimizer](https://www.gurobi.com/) in MADRID. |
|     `--volume=$HOME/madrid_local`      | Recommended |            The path to store local MADRID files at. This is not required, but recommended to negate data loss            |
|                `--name`                |     No      |                                What the container should be named for easy identification                                |
|                 `--rm`                 |     No      |                                   Remove the container after you are finished with it                                    |
|                 `-it`                  |     No      |                                             Run the container interactively                                              |
| `ghcr.io/helikarlab/madrid:CHOSEN-TAG` |     Yes     |                                           The container image and tag to start                                           |

## Accessing MADRID
Navigate to [http://localhost:8888](http://localhost:8888) in your browser

## Working in the Notebook
Open the Jupyter Notebook file, found at `/main/py/pipeline_paper_demo.ipynb` in the web interface.

Configuration files should be uploaded to `{{ site.data.terms.data_dir }}/config_files` and data files to `{{ site.data.terms.data_dir }}/data_matrices`, according to instructions in the notebook and provided templates

Update the file names in the Jupyter Notebook accordingly
