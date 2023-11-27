---
title: Knockout Simulation
sidebar: como_sidebar
permalink: como_args_knockout_simulation.html
folder: como
summary: "The arguments used in the `/main/src/knock_out_simulation.py` file"
last_updated: September 21, 2022
---

This file contains details about the Knockout Simulation arguments.

The default arguments for this script are as follows:
```python
cmd = ' '.join(['python3' , 'knock_out_simulation.py',
              '--context-model', '"{}"'.format(tissueSpecificModelfile),
              '--context-name', '"{}"'.format(context),
              '--disease-name', '"{}"'.format(dis),
              '--disease-up', '"{}"'.format(Disease_Up),
              '--disease-down', '"{}"'.format(Disease_Down),
              '--raw-drug-file', '"{}"'.format(drug_raw_file),
              '--reference-flux-file', '"{}"'.format(ref_flux_file),
              #'--test-all'
               ])
!{cmd}
```


{% include argument_table.html %}
