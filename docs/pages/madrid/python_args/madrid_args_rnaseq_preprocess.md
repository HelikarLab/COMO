---
title: RNASeq Preprocess Arguments
sidebar: madrid_sidebar
permalink: madrid_args_rnaseq_preprocess.html
folder: madrid
summary: "The arguments used in the `/main/py/rnaseq_preprocess.py` file"
last_updated: September 16, 2022
---

The default command for this block is as follows: 
```python
cmd = ' '.join(['python3', 'rnaseq_preprocess.py',
                '--context-names', '"{}"'.format(context_names),
                '--gene-format', '"{}"'.format(gene_format),
                '--taxon-id', '"{}"'.format(taxon_id),
                '--{}'.format(preprocess_mode)])
```

|                   Argument                    | Required? |                              Description                              |        Default        |
|:---------------------------------------------:|:---------:|:---------------------------------------------------------------------:|:---------------------:|
|               `--context-names`               |    Yes    |                         A list of cell types                          | `"['naiveB', 'smB']"` |
|                `--gene-format`                |    Yes    |                   The format your gene input is in                    |      `"Ensembl"`      |
|                 `--taxon-id`                  |    Yes    |                         The BioDBNet Taxon ID                         |       `"human"`       |
| `--create-matrix`<br>OR<br>`--provide-matrix` |    Yes    | Should the gene matrix be created, <br>or are you providing a matrix? |   `"create-matrix"`   |

