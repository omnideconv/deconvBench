#!/bin/bash
# With Java 21 on this cluster, copying a file larger than about 2 MB to an NFS share never returns,
# which blocks Nextflow when it publishes results. Java 17 is not affected.
export NXF_JAVA_HOME=${NXF_JAVA_HOME:-/usr/lib/jvm/java-17-openjdk-amd64}

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/configs/nextflow_norm_resolution.config run \
 /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf \
 -profile rectangle \
 -work-dir /nfs/scratch/nf-core_work/omnideconv.benchmark \
 -entry simulation_resolution_analysis "$@" \
 -resume