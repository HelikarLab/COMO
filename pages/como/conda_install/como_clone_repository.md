---
title: Clone the COMO Repository
permalink: como_clone_repository.html
summary: This guide will cover how to clone the COMO repository
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Clone the Repository

COMO is hosted at [github.com/HelikarLab/COMO](github.com/HelikarLab/COMO), and there are multiple methods of cloning.

## Using the Web

If you'd prefer to use a web browser, simply navigate to the [COMO repository](github.com/HelikarLab/COMO) and click the
green "Code" button. Then, click "Download ZIP" to download the repository as a ZIP file. Extract the ZIP file to a
location of your choice, and you've successfully cloned the repository

## Using Git

If you'd prefer to use Git, you can clone the repository using the following command:

{% include tip.html content="If you don't have Git installed, you can download
it <a href='https://git-scm.com/downloads' target='_blank'>here</a>" %}

```bash
git clone github.com/HelikarLab/COMO.git
```

## Using GitHub CLI

If you'd prefer to use the GitHub CLI, you can clone the repository using the following command:

{% include tip.html content="If you don't have the GitHub CLI installed, you can download
it <a href='https://cli.github.com/' target='_blank'>here</a>" %}

```bash
gh repo clone HelikarLab/COMO
```

## Summary

Now that you've cloned the repository, you're ready to install COMO's dependencies.
