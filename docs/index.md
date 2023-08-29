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

## Choosing an installation method

We provide two methods of installation: Conda/Mamba and Docker. Both methods are equally valid, but Docker can be more
difficult to set up due to the nature of mapping paths from the docker container to the host. If you are new to Docker,
we recommend using Conda/Mamba.

- [Install using Conda/Mamba](/como_conda_overview.html)
- [Install using Docker](/como_docker_overview.html)

## Accessing COMO

Navigate to [http://localhost:8888](http://localhost:8888) in your browser

## Working in the Notebook

Open the Jupyter Notebook file, found at `COMO.ipynb` in the web interface.

Configuration files should be uploaded to `{{ site.data.terms.data_dir }}/config_files` and data files
to `{{ site.data.terms.data_dir }}/data_matrices`, according to instructions in the notebook and provided templates

Update the file names in the Jupyter Notebook accordingly
