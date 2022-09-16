---
title: Arguments Overview
permalink: madrid_arguments_overview.html
tags: [getting_started]
summary: This is an overview of how arguments work in MADRID
sidebar: madrid_sidebar
last_updated: September 15, 2022
---

## What are Arguments?
Arguments in MADRID are the primary method of modifying settings. They are identical to command line arguments. The first argument is always the name of the script to be run. The rest of the arguments are the arguments that you want to pass to your MADRID program.

For example, if we wanted to enter the `--taxon-id` into a script, we would use the following format:

```bash
python3 script_name.py --taxon-id 9606
```

|     Argument     |             Description              |
|:----------------:|:------------------------------------:|
|     python3      | Execute python from the command line |
| `script_name.py` |  The name of the script to execute   |
| --taxon-id 9606  |  The taxon ID to be used in MADRID   |


Arguments of each script are listed in the sidebar, under the appropriate script name.  
[Python Arguments](madrid_args_python.html)  
[R Arguments](madrid_args_r.html)
