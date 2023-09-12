---
title: Conda Installation Overview
permalink: como_conda_installation.html
summary: This is how to install COMO without Docker
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Install Conda

The following instructions are taken directly from the Conda documentation. Please visit
the [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation)
for more details

### Windows

These three commands quickly and quietly install the latest 64-bit version of the installer and then clean up after
themselves. To install a different version or architecture of Miniconda for Windows, change the name of the .exe
installer in the curl command.

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe
start /wait "" miniconda.exe /S
del miniconda.exe
```

After installing, open the Conda prompot (miniconda3) program to use Miniconda3. The following commands initialize conda
for the `cmd.exe` and powershell shells:

```bash
conda init cmd.exe
conda init powershell
```

### MacOS

These four commands quickly and quietly install the latest M1 macOS version of the installer and then clean up after
themselves. To install a different version or architecture of Miniconda for macOS, change the name of the .sh installer
in the curl command.

```bash
mkdir -p ~/miniconda3
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:

```bash
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

### Linux

These four commands quickly and quietly install the latest 64-bit version of the installer and then clean up after
themselves. To install a different version or architecture of Miniconda for Linux, change the name of the .sh installer
in the wget command.

```bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

After installing, initialize your newly-installed Miniconda. The following commands initialize for bash and zsh shells:

```bash
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

## Initializing Conda

Once Conda is installed, you will need to initialize it. Initialization is important because it allows you to use the
`conda` command in your terminal. Initialization is done by running the following command:

```bash
conda init $SHELL
```

## Summary

Now that conda is installed, we can move to the next step of [cloning the COMO repository](/como_clone_repository.html)

