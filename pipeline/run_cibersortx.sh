#!/bin/bash

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/nextflow_norm.config run /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf -profile cibersortx -work-dir /localscratch/omnideconv_work -resume