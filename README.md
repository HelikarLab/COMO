# COMO: Constraint-based Optimizations of Metabolic Objectives

[![Unit Tests](https://github.com/HelikarLab/COMO/actions/workflows/unit_tests.yml/badge.svg)](https://github.com/HelikarLab/COMO/actions/workflows/unit_tests.yml)
[![Documentation](https://github.com/HelikarLab/COMO/actions/workflows/pages/pages-build-deployment/badge.svg?branch=master)](https://github.com/HelikarLab/COMO/actions/workflows/pages/pages-build-deployment)

This is the home page for the COMO pipeline.

[Create and view Issues](https://github.com/HelikarLab/COMO/issues)

## Quick Start

- [Install Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - Preferably, [install mamba](https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install) instead.
      Mamba is much faster than Conda and offers the same features
- [Clone this repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
    - `git clone https://github.com/HelikarLab/COMO.git`
- Change directories into the newly cloned repository
    - `cd COMO`
- Create a new conda environment
    - `conda env create -f environment.yaml`, <ins>**OR**</ins>
    - `mamba env create -f environment.yaml`
- Activate the new environment
    - `conda activate como`, <ins>**OR**</ins>
    - `mamba activate como`
- **IMPORTANT**: Install our modified version of zFPKM to allow for filtering insignificant local maxima during RNA-seq
  processing
    - `R -e "devtools::install_github('babessell1/zFPKM')"`
- Start the notebook server
    - `cd main && jupyter notebook` (for "retro" jupyter notebook look and feel), <ins>**OR**</ins>
    - `cd main && jupyter lab` (for the newer jupyter lab look and feel)

This will open a web browser with the Jupyter Notebook/Lab interface. From here, you can open the `COMO.ipynb` notebook
to get started

## Detailed Start

For more detailed information, please [view our documentation](https://helikarlab.github.io/COMO)

## Flow Charts

Please [follow this link](https://helikarlab.github.io/COMO/como_flowcharts.html) for flow charts

## Resources

Resources for packages used in COMO and other useful links, please
see [here](https://helikarlab.github.io/COMO/como_resources.html)
