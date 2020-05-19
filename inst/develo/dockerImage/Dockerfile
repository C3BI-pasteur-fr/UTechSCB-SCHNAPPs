# https://hub.docker.com/r/rocker/shiny
# https://github.com/rocker-org/shiny

### usage:
# docker run --rm -p 3838:3838 rocker/shiny
# docker run --rm -p 3838:3838 rocker/shiny:3.4.4
# docker run --user shiny --rm -p 3838:3838 schnapps_epdc:1.3 


# load R image
FROM r-base:3.6.3

# set working directory
WORKDIR /root

# install libraries needed
RUN apt-get -y update
RUN apt-get install -y \
    libxml2-dev \
    libgit2-dev  \
    libssl-dev  \
    libssh2-1-dev \
    libpoppler-cpp-dev \
    orca\
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "update.packages()" && \
    R -e "install.packages('devtools')" && \
    R -e "devtools::install_github('mul118/shinyMCE')" && \
    R -e "install.packages(c('BiocManager', 'pdftools'))" && \
    R -e "BiocManager::install(c('BiocSingular', 'SingleR'))" && \
    R -e "devtools::install_github('C3BI-pasteur-fr/UTechSCB-SCHNAPPs', dependencies = TRUE)" && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN', dependencies = TRUE)" && \
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

