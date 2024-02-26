FROM ubuntu:20.04

LABEL description="Image for Omnideconv benchmarking pipeline"
LABEL maintainer="Alexander Dietrich"

# needed so that R r-base installation does not get stuck
ENV DEBIAN_FRONTEND noninteractive

# install system dependencies
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y git \
    wget \
    vim \
    software-properties-common \
    dirmngr \
    gdebi \
    curl

# install R 4.2.3 https://docs.posit.co/resources/install-r/
RUN curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-4.2.3_1_amd64.deb
RUN gdebi r-4.2.3_1_amd64.deb -n
RUN ln -s /opt/R/4.2.3/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.2.3/bin/Rscript /usr/local/bin/Rscript

# install miniconda into /root/miniconda3
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /root/.local/share/r-miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/root/.local/share/r-miniconda/bin:${PATH}"
ARG PATH="/root/.local/share/r-miniconda/bin:${PATH}"

RUN conda --version


# create new conda environment, in which scaden and autogenes is installed
RUN conda create -n r-omnideconv python=3.10
RUN echo "source activate r-omnideconv" > ~/.bashrc

ENV PATH="/root/.local/share/r-miniconda/envs/omnideconv/bin:${PATH}"
ARG PATH="/root/.local/share/r-miniconda/envs/omnideconv/bin:${PATH}"

RUN conda run --no-capture-output -n r-omnideconv pip install scaden
RUN conda run --no-capture-output -n r-omnideconv pip install autogenes
RUN conda run --no-capture-output -n r-omnideconv pip install loguru
RUN conda run --no-capture-output -n r-omnideconv pip install pydeseq2==0.4.1  
RUN conda run --no-capture-output -n r-omnideconv pip install quadprog==0.1.11

# install omnideconv with all dependencies
RUN apt-get update -y && apt-get --no-install-recommends --fix-broken install -y curl \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgdal-dev \
    cmake


RUN R -e "install.packages(c('pak','pkgdepends','devtools'), repos='https://cloud.r-project.org')"
RUN R -e "pak::pkg_install(c('remotes','harmony','terra'))"
RUN echo  "GITHUB_PAT=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" >> /root/.Renviron
RUN R -e "remotes::install_github('omnideconv/omnideconv@benchmark', dependencies = TRUE, auth_token='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')"
RUN R -e "remotes::install_version('NMF', '0.26', repos='https://cloud.r-project.org')"


# set systems variable RETICULATE_PYTHON to new conda env, so that omnideconv will use it
RUN echo "RETICULATE_PYTHON = '/root/.local/share/r-miniconda/envs/r-omnideconv/bin/python'" > /root/.Rprofile
RUN R -e "reticulate::use_miniconda(condaenv = 'r-omnideconv', required = TRUE)"

# install benchmarking related dependencies
RUN R -e "install.packages(c('conflicted','docopt','readxl','tidyverse','profvis'), repos='https://cloud.r-project.org')"
RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('SimBu')"

# install Docker, so that CibersortX can run
RUN mkdir -p /etc/apt/keyrings
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /etc/apt/keyrings/docker.gpg
RUN echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
RUN apt-get update -y
RUN apt-get install -y --no-install-recommends \
    docker-ce \
    docker-ce-cli \
    containerd.io
