---
title: Proteomics Gen
sidebar: como_sidebar
permalink: como_args_proteomics_gen.html
folder: como
summary: "The arguments used in the `/main/src/proteomics_gen.py` file"
last_updated: September 21, 2022
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

{% include argument_table.html %}
