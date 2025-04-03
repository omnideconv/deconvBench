FROM eddelbuettel/r2u:20.04

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

# install miniconda into /root/miniconda3
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /root/.local/share/r-miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/root/.local/share/r-miniconda/bin:${PATH}"
ARG PATH="/root/.local/share/r-miniconda/bin:${PATH}"

RUN conda --version

# create r-omnideconv environment based on yml file
COPY environment.yml .
RUN conda env create -f environment.yml
RUN echo "source activate r-omnideconv" > ~/.bashrc


ENV PATH="/root/.local/share/r-miniconda/envs/omnideconv/bin:${PATH}"
ARG PATH="/root/.local/share/r-miniconda/envs/omnideconv/bin:${PATH}"

# install omnideconv with all dependencies
RUN apt-get --no-install-recommends --fix-broken install -y libcurl4-openssl-dev \
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


RUN R -e "install.packages(c('pak','pkgdepends','devtools','remotes','harmony','terra','igraph', 'RSQLite'), repos='https://cloud.r-project.org')"
RUN echo "GITHUB_PAT=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" >> /root/.Renviron
RUN R -e "remotes::install_github('omnideconv/omnideconv@benchmark', dependencies = TRUE, auth_token='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' )"
#RUN R -e "pak::pkg_install('omnideconv/SCDC')"
#RUN R -e "pak::pkg_install('omnideconv/BayesPrism')"
#RUN R -e "remotes::install_version('NMF', '0.27', repos='https://cloud.r-project.org')"


# set systems variable RETICULATE_PYTHON to new conda env, so that omnideconv will use it
RUN echo "RETICULATE_PYTHON = '/root/.local/share/r-miniconda/envs/r-omnideconv/bin/python'" > /root/.Rprofile
RUN R -e "reticulate::use_miniconda(condaenv = 'r-omnideconv', required = TRUE)"

# install benchmarking related dependencies
RUN R -e "install.packages(c('conflicted','docopt','readxl','tidyverse','profvis'), repos='https://cloud.r-project.org')"
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


