#!/bin/bash

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/configs/nextflow_norm_resolution.config run \
 /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf \
 -profile rectangle \
 -work-dir /nfs/scratch/nf-core_work/omnideconv.benchmark \
 -entry simulation_resolution_analysis "$@" \
 -resume