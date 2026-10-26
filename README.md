# deconvBench

[![CI](https://github.com/omnideconv/deconvBench/actions/workflows/ci.yml/badge.svg)](https://github.com/omnideconv/deconvBench/actions/workflows/ci.yml)

A Nextflow pipeline to benchmark second-generation cell-type deconvolution methods, most of them available through the R package [omnideconv](https://github.com/omnideconv/omnideconv). It deconvolves real and simulated bulk RNA-seq data with single-cell references and measures how well each method recovers the true cell-type fractions. The repository contains everything needed to reproduce the results of our benchmarking study.

## Contents

| Folder | Contents |
|---|---|
| [`pipeline/`](pipeline/) | The Nextflow pipeline: eight workflows (deconvolution of real bulks, simulations, reference subsampling, missing cell types, unknown content, spillover, cell-type resolution, sequencing technology), configs and tests. See [pipeline/README.md](pipeline/README.md). |
| [`docker/`](docker/) | Dockerfile of the standard image, and scripts that derive the rectangle and SECRET images and publish the images. |
| [`data/`](data/) | Small example datasets: `sc1`, `bulk1`, `bulk2`, and the CI test data in `data/ci/`. |
| [`visualisation/`](visualisation/) | Scripts behind the figures of the study. See [visualisation/README.md](visualisation/README.md). |
| [`.github/workflows/`](.github/workflows/) | Continuous integration: every workflow and every method on small test data. |

## Methods

AutoGeneS, BayesPrism, Bisque, CIBERSORTx, DWLS, MuSiC, Scaden and SCDC through omnideconv, plus rectangle, the hierarchical MuSiC (`music-hc`) and SECRET. Each runs in one of four images:

| Image | Methods | Source |
|---|---|---|
| standard 2.0 | BayesPrism, Bisque, CIBERSORTx, DWLS, MuSiC, music-hc, Scaden, SCDC | `alexd13/omnideconv_benchmark:2.0` (Docker), `oras://ghcr.io/omnideconv/deconvbench-images:standard-2.0` (Apptainer) |
| autogenes 1.3 | AutoGeneS | `alexd13/omnideconv_benchmark:1.3`, `oras://ghcr.io/omnideconv/deconvbench-images:autogenes-1.3` |
| rectangle 2.5 | rectangle (rectanglepy 1.7.0) | `oras://ghcr.io/omnideconv/deconvbench-images:rectangle-2.5` |
| secret 2.4 | SECRET | not public; build it with `docker/add_secret.sh` |

CIBERSORTx needs a token from [cibersortx.stanford.edu](https://cibersortx.stanford.edu/) and runs in its own Docker image, which is not part of our images. Details are in [pipeline/README.md](pipeline/README.md#images).

## Quick start

Requirements: Linux, [Nextflow](https://www.nextflow.io/) 25.10 with Java 17, and [Apptainer](https://apptainer.org/).

Deconvolve the two test bulks in `data/ci/` with MuSiC:

```bash
git clone https://github.com/omnideconv/deconvBench.git && cd deconvBench
apptainer pull standard.sif oras://ghcr.io/omnideconv/deconvbench-images:standard-2.0
cp -r data/ci run_input
export IMG_STANDARD=$PWD/standard.sif CI_INPUT=$PWD/run_input CI_OUT=$PWD/run_results CI_BIND=$PWD
nextflow -C pipeline/tests/ci/ci.config run pipeline/main.nf
```

The estimated cell-type fractions are in `run_results/results_general/*/deconvolution.rds`, the metrics against the ground truth in `results_metric.rds`. To run your own data, write a config as described in [pipeline/README.md](pipeline/README.md#quickstart).

## License

GPL-3.0, see [LICENSE](LICENSE). The methods keep their own licences; the SECRET repository states none, which is why its image is not public.
