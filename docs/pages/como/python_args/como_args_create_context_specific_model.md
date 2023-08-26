---
title: Create Context Specific Model
sidebar: como_sidebar
permalink: como_args_create_context_specific_model.html
folder: como
summary: "The arguments used in the `/main/src/create_context_specific_model.py` file"
last_updated: September 21, 2022
---

This script creates a context specific model

The default arguments for this script are as follows:
```python
cmd = ' '.join(['python3', 'create_context_specific_model.py', 
                  '--context-name', '"{}"'.format(key),
                  '--reference-model-file', '"{}"'.format(general_model_filepath), 
                  '--active-genes-filepath', '"{}"'.format(active_genes_filepath),
                  '--objective', '"{}"'.format(objective),
                  '--boundary-reactions-filepath', '"{}"'.format(boundary_rxns_filepath),
                  #'--exclude-reactions-filepath', '"{}"'.format(exclude_rxns_file),
                  '--force-reactions-filepath', '"{}"'.format(force_rxns_filepath),
                  '--algorithm', '"{}"'.format(recon_algorithm),
                  '--low-threshold', '"{}"'.format(low_thresh),
                  '--high-threshold', '"{}"'.format(high_thresh),
                  '--solver', '"{}"'.format(solver),
                  '--output-filetypes', '"{}"'.format(output_filetypes)])

!{cmd}
```

## Arguments

{% include argument_table.html %}
