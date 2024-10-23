FROM jupyter/minimal-notebook:latest AS builder

# Install UV for project dependencies
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/
COPY --chown=1000:100 main "${HOME}/main"
COPY --chown=1000:100 pyproject.toml "${HOME}"

RUN uv sync && \
    jupyter lab --generate-config && \
    echo "" > "${HOME}/.jupyter/jupyter_lab_config.py" && \
    echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_lab_config.py" && \
    echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_lab_config.py" && \
    echo "c.ServerApp.allow_root = True" >> "${HOME}/.jupyter/jupyter_lab_config.py" && \
    echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_lab_config.py" && \
    echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_lab_config.py"


VOLUME "${HOME}/main/data/local_files"