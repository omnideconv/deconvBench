# deconvBench

This repository serves as the main resource to reproduce our results from the benchmarking study of second-generation deconvolution methods. 

The repository is divided into four sections:
1) `docker/`

    Here we store the Dockerfile that can be used in the benchmarking pipeline. It includes the [omnideconv](https://github.com/omnideconv/omnideconv) R package with all of its dependencies. An image created by this Dockerfile can also easily be accessed from [Dockerhub](https://hub.docker.com/repository/docker/alexd13/omnideconv_benchmark/general). 
2) `pipeline/`

    This folder includes the nextflow pipeline that can be used to reproduce all results from our study and can also be extended to benchmark new deconvolution methods. A detailed description is included in the directory itself.
3) `visualisation/`

    This folder includes the scripts that generate the figures of our study. They directly use the results of the benchmarking pipeline as input; `visualisation/rectangle_revision/` holds the figures of the rectangle revision.
4) `data/`

    A small example dataset (`sc1`, `bulk1`, `bulk2`) that the tests and the quickstart use, and `prepare_tcga_brca.R`, which builds the input of the TCGA-BRCA run.

 