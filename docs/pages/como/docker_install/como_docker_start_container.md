---
title: Starting the Docker Container
permalink: como_docker_start_container.html
summary: This guide will cover how to start the COMO container
sidebar: como_sidebar
last_updated: August 29, 2023
---

## Before Starting

Before starting the container, several key components must be covered

- Mapping Paths
- Mapping Ports

### Mapping Paths

An important part of using Docker in general is mapping paths and ports from the host machine to the container. This is
done using the `--mount` and `--port` flags respectively. The `--mount` flag is used to map paths from the host machine
to the
container and the `--port` flag is used to map ports from the host machine to the container.

Becuase docker containers are isolated from the host machine, you must map paths from the host machine to the container

For example, if we want to provide a directory on our computer to COMO inside the container that holds data for us to
use when building a model,we would use the `--mount` flag to "map" the path from our computer to inside container.

```bash
docker run --mount /path/to/host/directory:/home/jovyan/main/data/local_files my_container_name
```

### Mapping Ports

Similar to mapping paths, we must also map ports from the host machine to the container. This is done using the `--port`
flag

For example, if we want to access the Jupyter Notebook running inside the container, we would map the port that the
notebook is running on to a port on our host machine. This is done using the `--port` flag

```bash
docker run --port 8888:8888 my_container_name
```

### Combining Flags

You can combine the `--mount` and `--port` flags to map paths and ports at the same time

```bash
docker run --mount /path/to/host/directory:/home/jovyan/main/data/local_files --port 8888:8888 my_container_name
```

## Starting the Container

Once you have determined the paths and ports you want to use, you can start the container using the `docker run` command

We are going to define a couple variables - `GRB_LICENSE_FILE` and `LOCAL_FILES` - to make the final command easier to
read. The descriptions of each component in this command is described in the table below

```bash
TAG="latest"  # Change this to the version you want to use
GRB_LICENSE_FILE=/path/to/gurobi.lic
LOCAL_FILES=/path/to/host/directory

docker run \
  --mount type=bind,source="${GRB_LICENSE_FILE}",target=/home/jovyan/gurobi.lic,readonly \
  --volume="${LOCAL_FILES}":/home/jovyan/main/data/local_files \
  --port 8888:8888  \
  --name como \
  --detach \
  ghcr.io/helikarlab/como:"${TAG}"
```

| Option                    | Required? | Description                                                                                                                                                                              |
|---------------------------|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `TAG`                     | No        | This is the version of COMO you want to use. If you do not change the variable to a specific version, the latest version will be used                                                    |
| `GRB_LICENSE_FILE`        | Yes       | This is the path to your Gurobi license file. This is required to use Gurobi inside the container                                                                                        |
| `LOCAL_FILES`             | Yes       | This is the path to the directory on your host machine that you want to map to the container. This is where you will store your data for COMO to use                                     |
| `--mount`                 | No        | This is used to tell COMO where the gurobi license is at, which is only required if you would like to use Gurobi                                                                         |
| `--volume`                | No        | This is used to tell COMO where your "local files" are at that you would like to access inside the container. It is only required if you would like to save your work after running OCMO |
| `--port`                  | Yes       | This is required to be able to access COMO in a web browser                                                                                                                              |
| `--name`                  | No        | This simply provides a name to the contianer so you know what is running                                                                                                                 |
| `--detach`                | No        | This allows you to close the terminal after executing the `docker run...` command. If you do not include this, COMO will close after closing the terminal                                |
| `ghcr.io/helikarlab/como` | Yes       | This tells Docker what image you would like to execute                                                                                                                                   |























