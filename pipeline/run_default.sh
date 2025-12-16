#!/bin/bash

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/nextflow_norm.config run /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf \
 -profile apptainer \
 -work-dir /nfs/scratch/nf-core_work/omnideconv.benchmark \
 -with-trace "$@"