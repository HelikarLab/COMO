FROM jupyter/r-notebook:latest

ARG GRB_SHORT_VERSION=9.5
ARG GRB_VERSION=9.5.0

ENV HOME /home/jovyan

COPY build_scripts/mamba_install.txt "${HOME}"/
COPY build_scripts/pip_install.txt "${HOME}"/

# Give ownership to jovyan user
COPY --chown=1000:100 pipelines "${HOME}"/work

# Install python-related items
RUN pip install --requirement "${HOME}/pip_install.txt" \
    && mamba install --yes --channel bioconda --channel conda-forge --file "${HOME}/mamba_install.txt" \
    && mamba clean --all --force-pkgs-dirs --yes \
    && jupyter trust "${HOME}/work/py/pipeline.ipynb" \
    && rm -f "${HOME}/pip_install.txt" \
    && rm -f "${HOME}/mamba_install.txt"

# Install gurbori
RUN wget --quiet https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xf gurobi${GRB_VERSION}_linux64.tar.gz \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* gurobi \
    && rm -rf gurobi/linux64/docs \

# Update jupyter notebook configuration \
RUN echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.root_dir = '${HOME}/work'" >> "${HOME}/.jupyter/jupyter_notebook_config.py"
