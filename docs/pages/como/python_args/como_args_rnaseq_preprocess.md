---
title: RNASeq Preprocess Arguments
sidebar: como_sidebar
permalink: como_args_rnaseq_preprocess.html
folder: como
summary: "The arguments used in the `/main/src/rnaseq_preprocess.py` file"
last_updated: September 21, 2022
---

The default command for this block is as follows: 
```python
cmd = ' '.join(['python3', 'rnaseq_preprocess.py',
                '--context-names', '"{}"'.format(context_names),
                '--gene-format', '"{}"'.format(gene_format),
                '--taxon-id', '"{}"'.format(taxon_id),
                '--{}'.format(preprocess_mode)])
```

{% include argument_table.html %}
