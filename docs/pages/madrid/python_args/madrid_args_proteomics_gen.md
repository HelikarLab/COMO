---
title: Proteomics Gen
sidebar: madrid_sidebar
permalink: madrid_args_proteomics_gen.html
folder: madrid
summary: "The arguments used in the `/main/py/proteomics_gen.py` file"
last_updated: September 16, 2022
---

This file creates proteomics matrices for downloaded proteomics data.

<br>

The default command for this block is as follows:
```python
cmd = ' '.join(['python3', 'proteomics_gen.py', 
      '--config-file', '"{}"'.format(proteomics_config_file),
      '--replicate-ratio', '"{}"'.format(rep_ratio),
      '--high-replicate-ratio', '"{}"'.format(high_rep_ratio),
      '--batch-ratio', '"{}"'.format(batch_ratio),
      '--high-batch-ratio', '"{}"'.format(high_batch_ratio),
      '--quantile', '"{}"'.format(quantile)])
```

|         Argument         | Required? |                                                                     Description                                                                     |               Default               |
|:------------------------:|:---------:|:---------------------------------------------------------------------------------------------------------------------------------------------------:|:-----------------------------------:|
|     `--config-file`      |    Yes    |                                                         The configuration file name to use                                                          | `proteomics_data_inputs_paper.xlsx` |
|   `--replicate-ratio`    |    Yes    |                                 The ratio of replicates required for a gene to be considered active in that sample                                  |               `0.75`                |
| `--high-replicate-ratio` |    Yes    | The ratio of replicates to be considered "high confidence".<br>High confidence genes are considered expressed no matter what other data sources say |                `1.0`                |
|     `--batch-ratio`      |    Yes    |                                   The ratio of batches required for a gene to be considered active in that sample                                   |               `0.75`                |
|   `--high-batch-ratio`   |    Yes    |  The ratio of batches to be considered "high confidence".<br>High confidence genes are considered expressed no matter what other data sources say   |                `1.0`                |
|       `--quantile`       |    Yes    |                                                  The percentage to be considered for each quantile                                                  |                `25`                 |
