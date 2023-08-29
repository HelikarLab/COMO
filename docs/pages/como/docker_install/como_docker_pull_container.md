---
title: Downloading the Container
permalink: como_docker_pull_container.html
summary: This guide will cover how to select the correct tag and pull the container
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Choosing a tag

### What are tags?

Tags in Docker are a way to provide a "version" to a container. Each tag is a different version of the container. For
example, the `latest` tag is the most recent version of the container. The `latest` tag is the default tag, so i you do
not specify a tag, Docker will download the `latest` tag.

If you would like to install an older version of COMO, you can specify the tag for that version. For example, if you
would like to install version 1.0.0 of COMO, you can specify the `1.0.0` tag as
such: `docker pull ghcr.io/helikarlab/como:1.0.0`

### Downloading the `latest` tag

To download the most recent, stable version of COMO, simply run the following command on the terminal:

```bash
docker pull ghcr.io/helikarlab/como
```

### Summary

Now that the container is downloaded, you can continue
to [Mapping Paths](/como_docker_mapping_paths.html)
