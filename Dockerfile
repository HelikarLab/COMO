ARG UBUNTU_RELEASE=20.04
ARG RPY2_VERSION=master
#FROM gurobi/python:9.5.0
FROM rpy2/base-ubuntu:$RPY2_VERSION-$UBUNTU_RELEASE

ARG DEBIAN_FRONTEND=noninteractive
ENV CRAN_MIRROR=https://cloud.r-project.org \
    CRAN_MIRROR_TAG=-cran35
ENV JUPYTER_ENABLE_LAB=1
ARG TINI_VERSION=v0.19.0
ENV SHELL=/bin/bash \
    NB_USER=jupyteruser \
    NB_UID=1000
#force
USER root

RUN \
  apt-get update \
  && apt-get install -y libssl-dev \
  && apt-get install -y libxml2-dev

#COPY install_py_libs.sh /opt/install_py_libs.sh
#COPY install_r_libs.sh /opt/install_r_libs.sh
#COPY install_jupyter.sh /opt/install_jupyter.sh
#COPY setup_jupyter.sh /opt/setup_jupyter.sh

COPY install_py_libs.sh /opt/
COPY install_r_libs.sh /opt/
COPY install_jupyter.sh /opt/install_jupyter.sh
COPY setup_jupyter.sh /opt/setup_jupyter.sh
  
RUN \
  apt-get update -qq \
  && apt-get install -y curl \
  && apt-get remove -y --purge nodejs npm \
  && sh /opt/install_py_libs.sh \
  && curl -sL https://deb.nodesource.com/setup_14.x | sudo -E bash - \
  && apt-get install -y nodejs \
  && wget -O - https://dl.yarnpkg.com/debian/pubkey.gpg | apt-key add - \
  && apt-get update -qq \
  && apt-get install -y yarn \
  && npm install -g configurable-http-proxy \
  && rm -rf /var/lib/apt/lists/* \
  && useradd -m -s /bin/bash -N -u "${NB_UID}" "${NB_USER}" \
  && sh /opt/install_jupyter.sh \
  && echo "${NB_USER}" 'ALL=(ALL) NOPASSWD: /usr/bin/apt-get' >> /etc/sudoers \
  && wget --quiet https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini \
      -P /usr/local/bin/ \
  && chmod +x /usr/local/bin/tini \
  && sh /opt/setup_jupyter.sh \
  && echo "Add Jupyter scripts emerging as ad hoc interface" \
  && git clone --depth=1 https://github.com/jupyter/docker-stacks.git /tmp/docker-stacks \
  && cd /tmp/docker-stacks/base-notebook \
  && sed -e 's/jovyan/'"${NB_USER}"'/g' start.sh > /usr/local/bin/start.sh \
  && cp start-notebook.sh /usr/local/bin/ \
  && chmod +x /usr/local/bin/start-notebook.sh \
  && cp start-singleuser.sh /usr/local/bin/ \
  && mkdir -p /etc/jupyter/ \
  && cp jupyter_notebook_config.py /etc/jupyter/ \
  && rm -rf /tmp/docker-stacks
  
RUN sh /opt/install_r_libs.sh  
  
  
COPY pipelines/ /home/"${NB_USER}"/work/
RUN chown -R "${NB_USER}" /home/"${NB_USER}"/work/py/
RUN chown -R "${NB_USER}" /home/"${NB_USER}"/work/data/
RUN chown -R "${NB_USER}" /usr/local/lib/R/site-library
#RUN chmod 755 /usr/local/bin/start-notebook.sh


USER $NB_USER

ENV HOME /home/${NB_USER}
WORKDIR ${HOME}

EXPOSE 8888

#RUN chmod +x /usr/local/bin/start-notebook.sh
#ENTRYPOINT ["sh", "/usr/local/bin/start-notebook.sh", "--"]
ENTRYPOINT ["/usr/local/bin/tini", "--"]
CMD ["start-notebook.sh"]