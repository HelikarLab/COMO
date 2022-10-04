---
title: Disease Analysis
sidebar: madrid_sidebar
permalink: madrid_args_disease_analysis.html
folder: madrid
summary: "The arguments used in the `/main/py/disease_analysis.py` file"
last_updated: September 21, 2022
---

This file contains details about the Disease Analysis arguments.

The default arguments for this script are as follows:
```python
cmd = ' '.join(['python3', 'disease_analysis.py',
              '--context-name', '"{}"'.format(context_name),
              '--config-file', '"{}"'.format(disease_config_file),
              '--data-source', '"{}"'.format(data_source),
              '--taxon-id', '"{}"'.format(taxon_id)])
!{cmd}
```


{% include argument_table.html %}
