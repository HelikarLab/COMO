FROM jupyter/r-notebook:latest as builder

COPY /environment.yaml "${HOME}/"
COPY --chown=1000:100 main "${HOME}/main"

# Install python-related items
# Remove "python" from the pinned file so we can install our own version
RUN sed -i '/^python/d' /opt/conda/conda-meta/pinned && \
    echo "auto_activate_base: true" >> "${HOME}/.condarc" && \
    jupyter trust "${HOME}/main/COMO.ipynb" && \
    echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    conda config --quiet --add channels conda-forge && \
    conda config --quiet --add channels bioconda && \
    conda config --quiet --add channels r && \
    rm -rf "${HOME}/main/tests"  # Remove tests, they are not required for running COMO

# Update base environment
RUN mamba env update --name=base --file="${HOME}/environment.yaml" && \
    R -e "devtools::install_github('babessell1/zFPKM')"

FROM jupyter/r-notebook:latest

COPY --from=builder ${HOME} ${HOME}

RUN pip cache purge && \
    conda clean --all --yes --force-pkgs-dirs

EXPOSE 8888
VOLUME /home/joyvan/main/data/local_files
