FROM jupyter/r-notebook:latest

ARG GRB_SHORT_VERSION=10.0
ARG GRB_VERSION=10.0.0

# Set gurobi environment variables
ENV GUROBI_HOME "${HOME}/gurobi/linux64"
ENV PATH "$PATH:$GUROBI_HOME/bin"
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:$GUROBI_HOME/lib"

COPY /environment.yaml "${HOME}/environment.yaml"
COPY --chown=1000:100 main "${HOME}"/main

# Install cobamp, bioservices, and troppo manually, as problems arose when installing using pip
RUN git clone https://github.com/babessell1/cobamp.git \
    && git clone https://github.com/babessell1/troppo.git \
    && pip install ./cobamp \
    && echo "COBAMP DONE" \
    && pip install ./troppo \
    && echo "TROPPO DONE" \
    && rm -rf cobamp troppo \
    && pip cache purge

# Install python-related items
# R channel will be placed at the bottom
RUN conda config --quiet --add channels r \
    && conda config --quiet --add channels bioconda \
    && conda config --quiet --add channels conda-forge \
    && mamba env update --name=base --file="${HOME}/environment.yaml" \
    && R -e "devtools::install_github('babessell1/zFPKM')" \
    # Install gurbori
    && wget --quiet --directory-prefix="${HOME}" https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xf "${HOME}/gurobi${GRB_VERSION}_linux64.tar.gz" \
    && rm -f "${HOME}/gurobi${GRB_VERSION}_linux64.tar.gz" \
    && mv ${HOME}/gurobi* ${HOME}/gurobi \
    && rm -f "${HOME}/environment.yaml" \
    && rm -r "${HOME}/work" \
    && pip cache purge \
    && conda clean --all --yes --force-pkgs-dirs

# Update jupyter notebook configuration
RUN jupyter trust "${HOME}/main/COMO.ipynb" \
    && echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py"

VOLUME /home/joyvan/main/data/local_files
