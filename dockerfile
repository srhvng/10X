###################################################################################################################
# EA 10X S3 Docker                                                                                                #
#-----------------------------------------------------------------------------------------------------------------#
# File Name:                                                                                                      #
#   Dockerfile                                                                                                    #
# Programmer:                                                                                                     #
#   Sarah Vang (sarah.vang@q2labsolutions.com)                                                                    #
# Last Edit:                                                                                                      #
#   9/14/2020                                                                                                     #
# Server:                                                                                                         #
#   vsd_docker01                                                                                                  #
#-----------------------------------------------------------------------------------------------------------------#
# Description:                                                                                                    #
#    Dockerizing the 10X tools (seurat,scrublet, and singleR) to avoid the                                        #
# ADC environment issues. This dockerfile starts from a R-4.0.0 docker                                            #
# image, adds python, additional R packages, and establishes the main                                             #
# directory (./) as the main place for a user to run the program.                                                 #
#                                                                                                                 #
#sudo /usr/bin/docker build --rm --no-cache -t ea-pipeline.quintiles.net/tenx_sss_pipeline:v1.0.0 -f dockerfile . #
###################################################################################################################
 
FROM rstudio/r-base:xenial
ARG R_VERSION=4.0.0
ARG OS_IDENTIFIER=ubuntu-1604
 
# Install R
RUN wget https://cdn.rstudio.com/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/*
 
CMD ["R"]
 
RUN Rscript -e "sessionInfo()"
 
#--------------------------------------------------------------
# INSTALL PYTHON and PACKAGES
#--------------------------------------------------------------
RUN apt-get update \
    && apt-get -y install software-properties-common \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get -y install python3 \
    && apt -y install python3-pip
 
RUN pip3 install --upgrade pip && pip3 install wheel llvmlite==0.31.0 && pip3 install argparse pandas pyyaml scipy matplotlib numpy scikit-learn scikit-image numba cython annoy umap-learn scrublet --use-feature=2020-resolver
 
RUN echo "##### PYTHON VERSION #####"; python3 --version
 
#--------------------------------------------------------------
# INSTALL R and PACKAGES
#--------------------------------------------------------------
RUN apt-get update && \
apt-get install -y libcurl4-gnutls-dev libxml2-dev gfortran libssl-dev libpng12-dev libjpeg-dev zlib1g-dev; # Needed for proper compile
 
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "sessionInfo()"
RUN Rscript -e "install.packages('devtools',repos='http://cran.us.r-project.org')"
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('SingleR','org.Hs.eg.db','Biobase'))"
RUN Rscript -e "install.packages('getopt',repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('dplyr',repos='http://cran.us.r-project.org')"
RUN Rscript -e "packageurl <- 'http://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.12.8.tar.gz' ;install.packages(packageurl, repos=NULL, type='source')"
RUN Rscript -e "install.packages('Seurat')"
RUN Rscript -e "install.packages('hdf5r')"
RUN Rscript -e "install.packages('patchwork',repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('purrr',repos='http://cran.us.r-project.org')"
RUN Rscript -e "library(getopt);library(dplyr);library(SingleR);library(Seurat);library(patchwork)"
