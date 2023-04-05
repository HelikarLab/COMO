# COMO: Constraint-based Optimizations of Metabolic Objectives 

[![Continuous Delivery](https://github.com/HelikarLab/COMO/actions/workflows/continuous_delivery.yml/badge.svg)](https://github.com/HelikarLab/COMO/actions/workflows/continuous_delivery.yml)
[![Unit Tests](https://github.com/HelikarLab/COMO/actions/workflows/unit_tests.yml/badge.svg?branch=master&event=pull_request)](https://github.com/HelikarLab/COMO/actions/workflows/unit_tests.yml)
[![pages-build-deployment](https://github.com/HelikarLab/COMO/actions/workflows/pages/pages-build-deployment/badge.svg?branch=master)](https://github.com/HelikarLab/COMO/actions/workflows/pages/pages-build-deployment)

This is the home page for the COMO pipeline.

For more detailed information, please [view the documentation here](https://helikarlab.github.io/COMO)

[Create and view Issues](https://github.com/HelikarLab/COMO/issues)

## Quick Start
- [Install Docker](https://docs.docker.com/install/)
- `sudo docker login`
- `sudo docker pull ghcr.io/helikarlab/como:latest`
- [Obtain a gurobi license](https://www.gurobi.com/academia/academic-program-and-licenses/)


### Using `docker run`
```bash
# We are going to export the license file location to an environment variable for use when running the docker container
GRB_LICENSE_FILE=/your/gurobi/license/file/location/gurobi.lic
LOCAL_FILES=/path/to/your/local/files

sudo docker run \
  --cpus=6 \
  -p 8888:8888 \
  --mount type=bind,source="${GRB_LICENSE_FILE}",target=/home/jovyan/gurobi.lic,readonly \
  --volume="${LOCAL_FILES}":/home/jovyan/main/data/local_files \
  --name como \
  -it \
  "ghcr.io/helikarlab/como:latest"
```


### Using `docker compose`
- Install the [docker compose](https://docs.docker.com/compose/install/) plugin
- Create a `docker-compose.yml` file with the following contents:
```yaml

```

- Run the docker image
  - `sudo docker run -p 8888:8888 --volume=$HOME/gurobi.lic:/opt/gurobi/gurobi.lic:ro  -v /$HOME/LocalComo:/home/jovyan/work/data/local_files --rm --name como --rm -it ghcr.io/helikarlab/como:latest`
- Open [http://127.0.0.1:8888](http://127.0.0.1:8888) from your browser  
- In your jupyter notebook, open `COMO.ipynb`
- Upload your configuration files to `data/config_files` and data files to `data/data_matrices` according to the instructions in the notebook and provided templates
  - Update the file names in the jupyter notebook accordingly.
- Run the notebook step by step, or run the step(s) by your needs


## Flow Charts
Please [follow this link](https://helikarlab.github.io/COMO/como_flowcharts.html) for flow charts

## Resources
Resources for packages used in COMO and other useful links, please see [here](https://helikarlab.github.io/COMO/como_resources.html)
