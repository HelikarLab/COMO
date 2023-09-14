FROM jupyter/r-notebook:latest

COPY /environment.yaml "${HOME}/environment.yaml"
COPY --chown=1000:100 main "${HOME}"/main

# Install python-related items
# Remove "python" from the pinned file so we can install our own version
RUN sed -i '/^python/d' /opt/conda/conda-meta/pinned && \
    conda config --quiet --add channels conda-forge && \
    conda config --quiet --add channels bioconda && \
    conda config --quiet --add channels r && \
    # Update conda
    mamba env update --name=base --file="${HOME}/environment.yaml" && \
    R -e "devtools::install_github('babessell1/zFPKM')" && \
    # Trust the jupyter notebook
    jupyter trust "${HOME}/main/COMO.ipynb" && \
    echo "c.ServerApp.ip = '0.0.0.0'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.root_dir = '${HOME}/main'" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.token = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    echo "c.ServerApp.password = ''" >> "${HOME}/.jupyter/jupyter_notebook_config.py" && \
    # Purge cache information, reducing image size
    pip cache purge && \
    conda clean --all --yes --force-pkgs-dirs && \
    rm -f "${HOME}/environment.yaml"

VOLUME /home/joyvan/main/data/local_files
