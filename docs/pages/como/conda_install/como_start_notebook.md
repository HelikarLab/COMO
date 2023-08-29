---
title: Starting COMO
permalink: como_start_notebook.html
summary: This guide will cover how to start COMO
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Starting COMO

If your `como` environment from previous steps is not active, activate it now:

```bash
conda activate como
```

To start COMO, run the following command:

```bash
jupyter lab
```

{% include tip.html content="If you prefer the 'retro' look and feel of jupyter notebooks, use '`jupyter notebook`'
instead" %}

## Accessing COMO

In the terminal, you will see the jupyter notebook beginning to start. Towards the end, a message similar to the
following will appear

```bash
    To access the server, open this file in a browser:
        file:///Users/joshl/Library/Jupyter/runtime/jpserver-86937-open.html
    Or copy and paste one of these URLs:
        http://localhost:8888/lab?token=c41795408eba585be71e1a8b652c23b775a5aa51d8075283
        http://127.0.0.1:8888/lab?token=c41795408eba585be71e1a8b652c23b775a5aa51d8075283
```

Copy one of the URLs starting with `http://`, and paste it into your browser. From here, you are ready to start using
COMO!


