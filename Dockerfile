FROM jupyter/r-notebook:latest

ARG GRB_SHORT_VERSION=9.5
ARG GRB_VERSION=9.5.0
ARG PYTHON_VERSION=3.10
ARG HOME="/home/jovyan"
ARG ENVIRONMENT_FILE="${HOME}/environment.yaml"
ARG JUPYTER_NOTEBOOK="${HOME}/work/py/pipeline.ipynb"
ARG JUPYTER_CONFIG="${HOME}/.jupyter/jupyter_notebook_config.py"

COPY /build_scripts/pip_install.txt "${HOME}/pip_install.txt"
COPY /build_scripts/mamba_install.txt "${HOME}/mamba_install.txt"
COPY /environment.yaml "${ENVIRONMENT_FILE}"

# Give ownership to jovyan user
COPY --chown=1000:100 /work "${HOME}"/work

# Install python-related items
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --add channels r \
    # Remove python from pinned versions; this allows us to update python. From: https://stackoverflow.com/a/11245372/13885200 \
    && sed -i "s/^python 3.*//" /opt/conda/conda-meta/pinned \
    && mamba install --yes python=${PYTHON_VERSION} \
    && mamba install --file "${HOME}/mamba_install.txt" \
    && python3 -m pip install -r "${HOME}/pip_install.txt" \
    #&& mamba env update --name base --file "${ENVIRONMENT_FILE}" \
    && mamba clean --all --force-pkgs-dirs --yes \
    && R -e 'devtools::install_github("babessell1/zFPKM")' \
    # Install jupyter extensions
    && jupyter labextension install @jupyter-widgets/jupyterlab-manager \
    && jupyter labextension install escher \
    && jupyter trust "${JUPYTER_NOTEBOOK}" \
    && rm -f "${ENVIRONMENT_FILE}" "${HOME}/pip_install.txt" "${HOME}/mamba_install.txt"

# Install gurbori
RUN wget https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xf gurobi${GRB_VERSION}_linux64.tar.gz \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* gurobi \
    && rm -rf gurobi/linux64/docs

# Update jupyter notebook configuration \
RUN echo "c.ServerApp.ip = '0.0.0.0'" >> "${JUPYTER_CONFIG}" \
    && echo "c.ServerApp.root_dir = '${HOME}/work'" >> "${JUPYTER_CONFIG}" \
    && echo "c.ServerApp.token = ''" >> "${JUPYTER_CONFIG}" \
    && echo "c.ServerApp.password = ''" >> "${JUPYTER_CONFIG}"

VOLUME /home/joyvan/work
