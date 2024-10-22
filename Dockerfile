FROM jupyter/r-notebook:latest as builder

COPY environment.yaml "${HOME}/environment.yaml"

# Install base packages and libraries
RUN mamba config --quiet --add channels conda-forge && \
    mamba config --quiet --add channels bioconda && \
    mamba config --quiet --add channels r && \
    sed -i '/^python/d' /opt/conda/conda-meta/pinned && \
    echo "auto_activate_base: true" >> "${HOME}/.condarc" && \
    mamba install --file="${HOME}/environment.yaml" --yes && \
    R -e "devtools::install_github('babessell1/zFPKM')" && \
    pip cache purge


FROM jupyter/r-notebook as production

COPY --from=builder "${HOME}/environment.yaml" "${HOME}/environment.yaml"
COPY --chown=1000:100 main "${HOME}/main"

# Remove tests directory
RUN rm -rf "${HOME}/main/tests"

# Configure jupyter notebook server
ENV JUPYTER_TOKEN=""
ENV JUPYTER_PASSWORD=""
ENV JUPYTER_IP="0.0.0.0"
ENV JUPYTER_ROOT_DIR="${HOME}/main"

RUN jupyter trust "${HOME}/main/COMO.ipynb"

VOLUME /home/jovyan/main/data/local_files

USER 1000