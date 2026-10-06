#!/bin/bash

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/configs/nextflow_norm.config run \
 /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf \
 -profile rectangle \
 -work-dir /nfs/scratch/nf-core_work/omnideconv.benchmark \
 -with-trace "$@" 