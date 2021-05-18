# https://hub.docker.com/r/rocker/shiny
# https://github.com/rocker-org/shiny

### usage:
# cd ..../inst/develo/dockerImage
# docker build -t pf2dock/schnapps:0.2  . 
# docker push pf2dock/schnapps:0.1
# docker run --rm -p 3838:3838 rocker/shiny
# docker run --rm -p 3838:3838 rocker/shiny:3.4.4
# local directory mapped to /data
# docker run --user shiny --rm -p 3838:3838 -v ./:/data schnapps:1.6.28 


# load R image
FROM r-base:4.0.4

# set working directory
WORKDIR /root

# install libraries needed
#\
#libssh2-1-dev \

RUN apt-get update && apt-get -y --fix-missing install git 
RUN apt-get update && apt-get install -y \
    libxml2-dev
RUN apt-get update && apt-get install -y \
#    libgit2-dev  
RUN git --version
RUN apt-get -y update && apt-get -y install libssl-dev  \
    apt-utils \
    ssh-client 
RUN apt-get update && apt-get install -y \
    libpoppler-cpp-dev
RUN apt-get update && apt-get install -y \
    orca \
    pandoc \
    sudo \
    gdebi-core 
RUN apt-get update && apt-get install -y \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libxt-dev \
    xtail \
    libgsl-dev \
    libhdf5-serial-dev \
    wget

# # Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
     VERSION=$(cat version.txt)  && \
     wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
     gdebi -n ss-latest.deb && \
     rm -f version.txt ss-latest.deb && \
     . /etc/environment
RUN     R -e "update.packages()" && \
    R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('mul118/shinyMCE')" && \
    R -e "install.packages(c('BiocManager', 'pdftools', 'ggplotify'))" && \
    R -e "BiocManager::install(c('BiocSingular', 'SingleR'))" && \
    R -e "BiocManager::install(c('multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer'))" && \
    R -e "BiocManager::install(c('IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics'))" && \
    R -e "BiocManager::install(c('DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment', 'S4Vectors'))" && \
    R -e "devtools::install_version('spatstat', version = '1.64-1', repos = 'http://cran.us.r-project.org')" && \
    R -e "devtools::install_github(repo = 'satijalab/seurat', dependencies = TRUE, ref = 'develop', upgrade='never')" && \
    R -e "remotes::install_github('C3BI-pasteur-fr/TemporaFork')" && \
    R -e "devtools::install_github('C3BI-pasteur-fr/UTechSCB-SCHNAPPs', dependencies = TRUE)"
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN', dependencies = TRUE)" && \
    R -e "install.packages(c('shinycssloaders', 'network', 'igraph', 'mclust', 'shinyTree', 'shinydashboard', 'hdf5r'))" && \
    R -e "setRepositories(ind = c(1,2,3));devtools::install_github('RausellLab/CelliD', ref = 'legacy')" && \
    R -e "install.packages(c('forcats', 'kohonen', 'SCORPIUS', 'shinyBS', 'threejs', 'DT', 'shinythemes'))"
RUN R -e "devtools::install_github('Albluca/distutils') " && \
    R -e "devtools::install_github('Albluca/ElPiGraph.R')" && \
    R -e "BiocManager::install('BiocParallel')" && \
    R -e "remotes::install_github('briatte/ggnetwork')"
 
RUN git clone https://github.com/baj12/SCHNAPPsContributions.git && \
    rm -r SCHNAPPsContributions/Dummy/  && \
    mkdir /data && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server
    
EXPOSE 3838

RUN mkdir -p /var/log/shiny-server
RUN mkdir -p /srv/shiny-server/schnapps

COPY app.R /srv/shiny-server/schnapps/

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

COPY *.RData /srv/shiny-server/schnapps/

COPY intro.Rmd /srv/shiny-server/schnapps/

COPY shiny-server.sh /usr/bin/shiny-server.sh

RUN chown -R shiny.shiny /var/log/shiny-server
RUN chown -R shiny.shiny /srv/shiny-server/schnapps


CMD ["/usr/bin/shiny-server.sh"]

