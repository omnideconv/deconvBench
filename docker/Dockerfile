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
    dirmngr

# install R 4.2 https://www.r-bloggers.com/2022/08/installation-of-r-4-2-on-ubuntu-22-04-1-lts-and-tips-for-spatial-packages/ 
RUN apt-get purge r-base* r-recommended r-cran-*
RUN apt autoremove
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt install -y r-base r-base-core r-recommended r-base-dev

# install miniconda into /root/miniconda3
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /root/miniconda3 \
    && rm Miniconda3-latest-Linux-x86_64.sh 

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN conda --version


# create new conda environment, in which scaden and autogenes is installed
RUN conda create -n omnideconv python=3.9
RUN echo "source activate omnideconv" > ~/.bashrc

ENV PATH="/root/miniconda3/envs/omnideconv/bin:${PATH}"
ARG PATH="/root/miniconda3/envs/omnideconv/bin:${PATH}"

RUN pip install scaden
RUN pip install autogenes

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

RUN R -e "install.packages(c('pak','pkgdepends','devtools'))"
RUN R -e "install.packages(c('harmony','remotes'))"
RUN R -e "remotes::install_github('omnideconv/CDSeq')"
RUN R -e "pak::pkg_install('omnideconv/omnideconv', dependencies = TRUE)"
#RUN R -e "pak::pkg_install('omnideconv/omnideconv@scaden_additions')"


# set systems variable RETICULATE_PYTHON to new conda env, so that omnideconv will use it
RUN echo "RETICULATE_PYTHON = '/root/miniconda3/envs/omnideconv/bin/python'" > /root/.Rprofile

# install benchmarking related dependencies
RUN R -e "install.packages(c('conflicted','docopt','readxl','tidyverse'))"
RUN R -e "install.packages(c('BiocManager'))"
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

