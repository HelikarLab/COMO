FROM jupyter/r-notebook:latest

ARG COMO_VERSION
ENV COMO_VERSION=${COMO_VERSION}

# Set gurobi environment variables
ENV GUROBI_HOME "${HOME}/gurobi/linux64"
ENV PATH "$PATH:$GUROBI_HOME/bin"
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:$GUROBI_HOME/lib"

COPY /environment.yaml "${HOME}/environment.yaml"
COPY --chown=1000:100 main "${HOME}"/main

# Install python-related items
# Remove python from pinned versions; this allows us to update python. From: https://stackoverflow.com/a/11245372
RUN sed -i "s|^python .*||" /opt/conda/conda-meta/pinned \
    && mamba env update --quiet --name=base --file="${HOME}/environment.yaml" \
    && mamba clean --quiet --all --force-pkgs-dirs --yes \
    && R -e "devtools::install_github('babessell1/zFPKM')" \
    && rm -rf "${HOME}/environment.yaml"

# Update jupyter notebook configuration \
RUN jupyter trust "${HOME}/main/COMO.ipynb" \
    && echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" \
    && echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py"

VOLUME /home/joyvan/main
