# to start the latest:
# docker pull pf2dock/schnapps
# docker run --user shiny --rm -p 3838:3838 -v /Users/bernd/Rstudio/:/data pf2dock/schnapps
# in e.g. chrom:
# http://localhost:3838/

# https://hub.docker.com/r/rocker/shiny
# https://github.com/rocker-org/shiny

### usage:
# cd inst/develo/dockerImage
# docker build -t pf2dock/schnapps:0.2  . 
# docker push pf2dock/schnapps:0.1
# docker run --rm -p 3838:3838 rocker/shiny
# docker run --rm -p 3838:3838 rocker/shiny:3.4.4
# local directory mapped to /data
# docker run --user shiny --rm -p 3838:3838 -v ./:/data schnapps:1.6.28 
# docker run --user shiny --rm -p 3838:3838 -v /Users/bernd/Rstudio/:/data pf2dock/schnapps:0.2

# docker build -f DockerfileUpdate.txt -t pf2dock/schnapps:0.31 .
# docker run --user shiny --rm -p 3838:3838 -v /Users/bernd/Rstudio/:/data pf2dock/schnapps:0.31

# docker run --user shiny --rm -p 3838:3838 -v /Users/bernd/Rstudio/:/data registry-gitlab.pasteur.fr/bernd/utechscb-schnapps

# load R image
FROM r-base:4.3.0

# set working directory
WORKDIR /root

# install libraries needed
#\
#libssh2-1-dev \
RUN apt-get update && apt-get upgrade -y && apt-get -y --fix-missing install git && apt-get install  -y libxml2-dev \
&& apt-get install -y   libgit2-dev  
RUN git --version
RUN apt-get -y update && apt-get -y install libssl-dev \
   && apt-get install -y apt-utils \
   && apt-get install -y ssh-client \
   && apt-get install -y libpoppler-cpp-dev \
   && apt-get install -y  orca \
   && apt-get install -y    pandoc \
   && apt-get install -y    sudo \
   && apt-get install -y    gdebi-core 
RUN apt-get update && apt-get install -y     libcurl4-gnutls-dev \
   && apt-get install -y     libxt-dev \
   && apt-get install -y     xtail \
   && apt-get install -y     libxt-dev \
   && apt-get install -y     libgsl-dev \
   && apt-get install -y     libhdf5-serial-dev \
   && apt-get install -y     wget \
   && apt-get install -y     libcairo2-dev \
   && apt-get install -y     libcurl4-openssl-dev\
   && apt-get install -y procps\
   && apt install -y libharfbuzz-dev libfribidi-dev\
   && apt-get install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN sudo apt install cmake  cmake-data libarchive13 libjsoncpp25 librhash0 libuv1
# # Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
     VERSION=$(cat version.txt)  && \
     wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
     gdebi -n ss-latest.deb && \
     rm -f version.txt ss-latest.deb && \
     . /etc/environment
RUN R -e "update.packages()" && \
    R -e "install.packages('devtools') "
#
RUN R -e "install.packages(c('BiocManager', 'pdftools', 'ggplotify'))" && \
    R -e "BiocManager::install(c('BiocSingular', 'SingleR', 'GSVA', 'GSEABase'))" && \
    R -e "BiocManager::install(c('multtest', 'limma', 'Biobase', 'monocle', 'rtracklayer'))" && \
    R -e "BiocManager::install(c('IRanges', 'GenomeInfoDb', 'GenomicRanges', 'BiocGenerics'))" && \
    R -e "BiocManager::install(c('DESeq2', 'MAST', 'SingleCellExperiment', 'SummarizedExperiment', 'S4Vectors'))"  && \
    R -e "devtools::install_version('spatstat', version = '1.64-1', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_github(repo = 'satijalab/seurat', dependencies = TRUE)" 
RUN R -e "remotes::install_github('C3BI-pasteur-fr/TemporaFork')"
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN', dependencies = TRUE)" && \
     R -e "install.packages(c('shinycssloaders', 'network', 'igraph', 'mclust', 'shinyTree', 'shinydashboard', 'hdf5r'))"
RUN R -e "install.packages(c('forcats', 'kohonen', 'SCORPIUS', 'shinyBS', 'threejs', 'DT', 'shinythemes'))"
RUN R -e "devtools::install_github('Albluca/distutils') " && \
     R -e "devtools::install_github('Albluca/ElPiGraph.R')" && \
     R -e "devtools::install_github('haowulab/Wind', build_opts = c('--no-resave-data'))" && \
     R -e "BiocManager::install('BiocParallel')" && \
     R -e "remotes::install_github('briatte/ggnetwork')" \
     R -e "devtools::install_github('theislab/destiny', dependencies = TRUE)"
RUN apt-get install lsof && \
    apt-get install net-tools && \
    apt-get update -y; sudo apt-get install -y r-cran-rcppeigen
RUN R -e "devtools::install_github('mul118/shinyMCE')"  && \
    R -e "devtools::install_github('haowulab/Wind', build_opts = c('--no-resave-data'))" && \
    R -e "install.packages(c('psychTools', 'pryr', 'colourpicker', 'Seurat', 'scater', 'shinyjs', 'tidySingleCellExperiment', 'shinyWidgets', 'dendsort', 'reactlog', 'shinydashboardPlus'))" && \
    R -e "BiocManager::install(c('scater', 'tidySingleCellExperiment'))" && \
    R -e "BiocManager::install(c('InteractiveComplexHeatmap'))"  && \
    R -e "install.packages(c('shinyjqui'))" && \
    R -e "install.packages(c('rintrojs', 'ggalluvial'))" && \
    R -e "devtools::install_github('nghiavtr/BPSC')" && \
    R -e "BiocManager::install('DEsingle')" && \
    R -e "devtools::install_github('nghiavtr/BPSC')" && \
    R -e "BiocManager::install('DESeq2')" && \
    R -e "BiocManager::install('edgeR')" && \
    R -e "BiocManager::install('MAST')" && \
    R -e "BiocManager::install('monocle')" && \
    R -e "BiocManager::install('limma')" && \
    R -e "BiocManager::install('Seurat')" && \
    R -e "devtools::install_github('statOmics/zingeR')" && \
    R -e "BiocManager::install('SingleCellExperiment')" && \
    R -e "BiocManager::install('scater')" && \
    R -e "devtools::install_github('Zhangxf-ccnu/scDEA')" && \
    R -e "install.packages('manhattanly')" && \
    R -e "install.packages('shinyAce')" && \
    R -e "install.packages('Hmisc')" && \
    R -e "install.packages('orca')" && \
    R -e "install.packages('debugme')" && \
    R -e "BiocManager::install('TRONCO')" && \
    R -e "BiocManager::install('SIMLR')"


RUN R -e " devtools::install_github('C3BI-pasteur-fr/UTechSCB-SCHNAPPs', dependencies = F )"


 #RUN git clone https://github.com/baj12/SCHNAPPsContributions.git && \
#     rm -r SCHNAPPsContributions/Dummy/  && \
#     mkdir /data && \
RUN    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
     chown shiny:shiny /var/lib/shiny-server
    
EXPOSE 3838

RUN mkdir -p /var/log/shiny-server
RUN mkdir -p /srv/shiny-server/schnapps

COPY app.R /srv/shiny-server/schnapps/

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

COPY *.RData /srv/shiny-server/schnapps/

COPY intro.Rmd /srv/shiny-server/schnapps/

COPY shiny-server.sh /usr/bin/shiny-server.sh
ENV SHINY_LOG_STDERR=1
RUN chown -R shiny:shiny /var/log/shiny-server
RUN chown -R shiny:shiny /srv/shiny-server/schnapps
RUN chown -R shiny:shiny /var/lib/shiny-server

CMD ["/usr/bin/shiny-server.sh"]

# CMD bash
