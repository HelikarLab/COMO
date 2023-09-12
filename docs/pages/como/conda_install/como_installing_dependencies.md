---
title: Installing Dependencies
permalink: como_install_dependencies.html
summary: This guide will cover how to install COMO's dependencies in a Conda environment
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Overview

We will be installing COMO's dependencies in a Conda (or Mamba!) virtual environment. If you don't know what a virtual
environment is, here is a brief explanation:

> A virtual environment is a self-contained and isolated workspace that allows you to manage and segregate dependencies
> for different projects. This means that you can have different sets of packages and libraries installed in separate
> virtual environments, preventing conflicts between projects that might require different versions of the same package.
> Virtual environments help keep your project dependencies organized, making it easier to maintain, share, and reproduce
> your software environments. When you work within a virtual environment, any packages you install or modify only affect
> that specific environment, leaving your system-wide configurations unchanged. This is especially useful when working
> on
> projects with specific version requirements or when collaborating with others, as it ensures consistent and controlled
> development environments.

## Creating a virtual environment

The first step is to create the actual virtual environment. This command should done on the command line

```bash
conda create -n como
```

## Activating the Environment

The environment has been created, but it hasn't been activated yet. Activate the environment by executing the below line

```bash
conda activate como
```

## Install Mamba

Mamba is much faster than Conda, due to it being written in C++. It is a drop-in replacement for Conda, so you can use
it in place of Conda without any issues. To install Mamba, execute the following command:

```bash
conda install --channel conda-forge mamba
```

## Installing Dependencies

Once the environment has been created, we can install the dependencies required by COMO.

{% include warning.html content="You must execute this in the root directory of the COMO folder that was downloaded in
the previous step" %}

```bash
# Change directories into the COMO folder
cd COMO/
mamba env update --file environment.yml
```

Once this is done, the final step is installing our customized version of zFPKM that allows for filtering insignificant
local maxima during RNA-seq processing. During the "Installing Dependencies" step, R was installed. We can use that now
to install zFPKM from our source

{% include important.html content="Installing zFPKM from our source is an important step in the installation process.
Failing to do so will result in errors when running COMO." %}

```bash
R -e "devtools::install_github('babessell1/zFPKM')"
```

## Summary

Now that you've installed COMO's dependencies, you're ready to [start using COMO](/como_start_notebook.html)!
