mkdir /home/jupyteruser/rlibs
#R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
#R -e 'BiocManager::install("affy", dep=TRUE, ask=FALSE)'
#R -e 'BiocManager::install("agilp", dep=TRUE, ask=FALSE)'
#R -e 'BiocManager::install("limma", dep=TRUE, ask=FALSE)'
#R -e 'BiocManager::install("hgu133acdf", dep=TRUE, ask=FALSE)'
#R -e 'install.packages(c("rzmq","repr","IRkernel","IRdisplay"), repos = c("http://irkernel.github.io/", getOption("repos")), type = "source")'
R -e 'install.packages(c("tidyverse", "sjmisc"), dependencies=TRUE, repos="'"${CRAN_MIRROR}"'", lib="/home/jupyteruser/rlibs")'
