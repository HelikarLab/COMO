---
title: RNASeq Preprocess Arguments
sidebar: madrid_sidebar
permalink: madrid_args_rnaseq_preprocess.html
folder: madrid
summary: "The arguments used in the `/main/py/rnaseq_preprocess.py` file"
last_updated: September 14, 2022
---

|                Argument                 | Required? |                              Description                              |
|:---------------------------------------:|:---------:|:---------------------------------------------------------------------:|
|            `--context-names`            |    Yes    |                         A list of cell types                          |
|             `--gene-format`             |    Yes    |                   The format your gene input is in                    |
|              `--taxon-id`               |    Yes    |                         The BioDBNet Taxon ID                         |
| `--create-matrix` OR `--provide-matrix` |    Yes    | Should the gene matrix be created, <br>or are you providing a matrix? |

