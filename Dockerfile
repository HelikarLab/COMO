FROM jupyter/r-notebook:latest

ARG GRB_SHORT_VERSION=10.0
ARG GRB_VERSION=10.0.0
ARG HOME="/home/joyvan"

# Set gurobi environment variables
ENV GUROBI_HOME "${HOME}/gurobi/linux64"
ENV PATH "$PATH:$GUROBI_HOME/bin"
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:$GUROBI_HOME/lib"

COPY /environment.yaml "${HOME}/environment.yaml"
COPY --chown=1000:100 main "${HOME}"/main

# Update jupyter notebook configuration
RUN jupyter trust "${HOME}/main/COMO.ipynb" \
    && echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py"

# Install python-related items
RUN conda config --quiet --add channels conda-forge \
    && conda config --quiet --add channels bioconda \
    && conda config --quiet --add channels r \
    && mamba env update --quiet --name=base --file="${HOME}/environment.yaml" \
    && mamba clean --quiet --all --force-pkgs-dirs --yes \
    && R -e "devtools::install_github('babessell1/zFPKM')" \
    # Install gurbori \
    && wget --quiet https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xf gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* "${HOME}/gurobi" \
    && rm -f "${HOME}/environment.yaml" \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz

VOLUME /home/joyvan/main/data/local_files
