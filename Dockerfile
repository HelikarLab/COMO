#### Base Image ####

ARG UBUNTU_RELEASE=20.04
FROM babessell/base-ubuntu:$RPY2_VERSION-$UBUNTU_RELEASE
#FROM babessell/rpy2-41:latest


#### Set Container Args ####

ARG DEBIAN_FRONTEND=noninteractive
ENV CRAN_MIRROR=https://cloud.r-project.org \
    CRAN_MIRROR_TAG=-cran41
	
ENV JUPYTER_ENABLE_LAB=1

ARG TINI_VERSION=v0.19.0

ENV SHELL=/bin/bash \
    NB_USER=jupyteruser \
    NB_UID=1000
	
USER root


#### Install Gurobi #####
# from https://github.com/Gurobi/docker-optimizer/blob/master/9.5.0/Dockerfile

ARG GRB_VERSION=9.5.0
ARG GRB_SHORT_VERSION=9.5
WORKDIR /opt

RUN apt-get update \
    && apt-get install --no-install-recommends -y\
       ca-certificates  \
       wget \
    && update-ca-certificates \
    && wget -v https://packages.gurobi.com/${GRB_SHORT_VERSION}/gurobi${GRB_VERSION}_linux64.tar.gz \
    && tar -xvf gurobi${GRB_VERSION}_linux64.tar.gz  \
    && rm -f gurobi${GRB_VERSION}_linux64.tar.gz \
    && mv -f gurobi* gurobi \
    && rm -rf gurobi/linux64/docs
	
	
#### Bioconductor dependencies ####

RUN \
  apt-get update \
  && apt-get install -y libssl-dev \
  && apt-get install -y libxml2-dev

#### Copy Install Scripts ####

COPY install_py_libs.sh /opt/
COPY install_r_libs.sh /opt/
COPY install_jupyter.sh /opt/install_jupyter.sh
COPY setup_jupyter.sh /opt/setup_jupyter.sh

#### Install Juypter and Python Libraries ####

# TODO: Split this into multiples parts that are more readable

RUN apt-get update -qq \
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
  && chmod +x /usr/local/bin/start.sh \
  && cp start-singleuser.sh /usr/local/bin/ \
  && mkdir -p /etc/jupyter/ \
  && cp jupyter_server_config.py /etc/jupyter/ \
  && rm -rf /tmp/docker-stacks
  
  
#### Gurobi Setup ####

LABEL vendor="Gurobi"
LABEL version=${GRB_VERSION}

RUN apt-get update \
    && apt-get install --no-install-recommends -y\
       ca-certificates  \
       p7zip-full \
       zip \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*
	
RUN python3 -m pip --no-cache-dir install gurobipy

WORKDIR /opt/gurobi
#COPY --from=buildoptimizer /opt/gurobi .

ENV GUROBI_HOME /opt/gurobi/linux64
ENV PATH $PATH:$GUROBI_HOME/bin
ENV LD_LIBRARY_PATH $GUROBI_HOME/lib
  
  
#### Install R Libraries ####
  
RUN sh /opt/install_r_libs.sh 


#### Ownership ####  

COPY pipelines/ /home/"${NB_USER}"/work/
RUN chown -R "${NB_USER}" /home/"${NB_USER}"/work/py/
RUN chown -R "${NB_USER}" /home/"${NB_USER}"/work/data/
RUN chown -R "${NB_USER}" /usr/local/lib/R/site-library
#RUN chmod 755 /usr/local/bin/start-notebook.sh


#### Set workspace and run Juypterlab ####

USER $NB_USER

ENV HOME /home/${NB_USER}
WORKDIR ${HOME}

EXPOSE 8888
ENTRYPOINT ["/usr/local/bin/tini", "--"]
CMD ["start-notebook.sh"] # rebuild
