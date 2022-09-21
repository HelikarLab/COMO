---
title: RNASeq Gen Arguments
sidebar: madrid_sidebar
permalink: madrid_args_rnaseq_gen.html
folder: madrid
summary: "The arguments used in the `/main/py/rnaseq_gen.py` file"
last_updated: September 21, 2022
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

{% include argument_table.html %}
