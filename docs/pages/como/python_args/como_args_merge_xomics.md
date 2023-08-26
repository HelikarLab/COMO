---
title: Merge Omics Gen
sidebar: como_sidebar
permalink: como_args_merge_xomics.html
folder: como
summary: "The arguments used in the `/main/src/merge_xomics.py` file"
last_updated: September 21, 2022
---

This file merges the omics data into a single matrix.

The default arguments for this script are as follows:
```python
cmd = ' '.join(['python3', 'merge_xomics.py', 
      '--merge-distribution',
      #'--microarray-config-file', '"{}"'.format(microarray_config_file),  # Only used if using microarray data
      '--total-rnaseq-config-file', '"{}"'.format(trnaseq_config_file),
      '--mrnaseq-config-file', '"{}"'.format(mrnaseq_config_file),
      #'--scrnaseq-config-file', '"{}"'.format(scrnaseq_config_file),  # Only used if using scRNA-seq data
      '--proteomics-config-file', '"{}"'.format(proteomics_config_file),
      '--expression-requirement', '"{}"'.format(expression_requirement),
      '--requirement-adjust', '"{}"'.format(requirement_adjust),
      '--total-rnaseq-weight', '"{}"'.format(tweight),
      '--mrnaseq-weight', '"{}"'.format(mweight),
      #'--single-cell-rnaseq-weight', '"{}"'.format(scweight),  # Only used if using scRNA-seq data
      '--protein-weight', '"{}"'.format(pweight),
      '--no-hc'
               ])
!{cmd}
```

## Arguments

{% include argument_table.html %}
