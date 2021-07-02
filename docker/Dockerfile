#Download base image ubuntu 18.04
FROM ubuntu:18.04
# Update Ubuntu Software repository
RUN apt-get update && \
        apt-get install -y software-properties-common && \
        apt-get update -y  && \
        apt-get install -y build-essential python3 python3-dev python3-pip apt-transport-https nano && \
        apt-get install -y git curl wget && \
        # update pip
        python3 -m pip install pip --upgrade && \
        python3 -m pip install wheel
WORKDIR /root
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt-get update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -y install tzdata
RUN apt-get install -y r-base
RUN apt-get install -y gdebi-core wget
RUN wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.5001-amd64.deb
RUN gdebi rstudio-server-1.2.5001-amd64.deb
RUN apt-get install -y python3-rpy2
RUN R -e "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")"
RUN R -e "BiocManager::install('affy', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('agilp', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('limma', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('hgu133acdf', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('edgeR', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('genefilter', dep=TRUE, ask=FALSE)"
RUN R -e "BiocManager::install('biomaRt', dep=TRUE, ask=FALSE)"
RUN R -e "install.packages(c('rzmq','repr','IRkernel','IRdisplay'), repos = c('http://irkernel.github.io/', getOption('repos')), type = 'source')"
RUN R -e "install.packages('sjmisc',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
RUN python3.6 -m pip install numpy --upgrade 
RUN R -e "IRkernel::installspec()"
COPY . .
#RUN /root/MCR_R2018a_glnxa64_installer/install -destinationFolder /root/mcr -agreeToLicense yes -mode silent
#RUN python3 /root/createTSM_GIMME/for_redistribution_files_only/setup.py install
#ENV XAPPLRESDIR=/root/mcr/v94/X11/app-defaults
#ENV LD_LIBRARY_PATH=/root/mcr/v94/runtime/glnxa64:/root/mcr/v94/bin/glnxa64:/root/mcr/v94/sys/os/glnxa64:/root/mcr/v94/extern/bin/glnxa64
#RUN rm -rf /root/MCR_R2018a_glnxa64_installer
RUN rm -rf rstudio-server-1.2.5001-amd64.deb
CMD ["rstudio-server", "start"]
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--allow-root"]
