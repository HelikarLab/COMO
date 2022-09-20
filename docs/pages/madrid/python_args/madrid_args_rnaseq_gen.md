---
title: RNASeq Gen Arguments
sidebar: madrid_sidebar
permalink: madrid_args_rnaseq_gen.html
folder: madrid
summary: "The arguments used in the `/main/py/rnaseq_gen.py` file"
last_updated: September 16, 2022
---

This file creates RNA sequencing matrices for total and mRNA preparation methods. As a result,this section is called twice, once with `--prep-method` set to `"total"`, and the second time set to `"mrna"`.

<br>

The default command for this block is as follows:
```python
cmd = ' '.join(['python3', 'rnaseq_gen.py', 
      '--config-file', '"{}"'.format(trnaseq_config_file), 
      '--replicate-ratio', '"{}"'.format(rep_ratio),   
      '--batch-ratio', '"{}"'.format(group_ratio),        
      '--high-replicate-ratio', '"{}"'.format(rep_ratio_h),    
      '--high-batch-ratio', '"{}"'.format(group_ratio_h),   
      '--filt-technique', '"{}"'.format(technique),  
      '--min-zfpkm', '"{}"'.format(min_zfpkm),
      '--library-prep', '"{}"'.format(prep_method)])       


```



|         Argument         | Required? |                                   Description                                    |               Default               |
|:------------------------:|:---------:|:--------------------------------------------------------------------------------:|:-----------------------------------:|
|     `--config-file`      |    Yes    |                        The configuration file name to use                        | `trnaseq_data_inputs_paper_rm.xlsx` |
|   `--replicate-ratio`    |    Yes    |    The proportion of replicates with expression required for high-confidence     |               `0.75`                |
|     `--batch-ratio`      |    Yes    |           The proportion of samples with expression required for gene            |               `0.75`                |
| `--high-replicate-ratio` |    Yes    |    The proportion of replicates with expression required for high-confidence     |                `1.0`                |
|   `--high-batch-ratio`   |    Yes    |    The proportion of replicates with expression required for high-confidence     |                `1.0`                |
|    `--filt-technique`    |    Yes    | The filtering technique; one of: (1) `"quantile"`, (2) `"cpm"`, or (3) `"zfpkm"` |              `"zfpkm"`              |
|      `--min-zfpkm`       |    Yes    |      The cutoff for CPM filtering, if using zFPKM as the `--filt-technique`      |                `-3`                 |
|     `--library-prep`     |    Yes    |  The library prepartion method; one of: (1) `"total"`, (2) `"mRNA"`, or `"SC"`   |              `"total"`              |
