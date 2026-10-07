#!/bin/bash
# With Java 21 on this cluster, copying a file larger than about 2 MB to an NFS share never returns,
# which blocks Nextflow when it publishes results. Java 17 is not affected.
export NXF_JAVA_HOME=${NXF_JAVA_HOME:-/usr/lib/jvm/java-17-openjdk-amd64}

# Rerun of the main workflow on the corrected (full NCBI) morandini bulk data with HaoCleaned-sampled.
# Usage: run_morandini_ncbi.sh <apptainer|autogenes|cibersortx> [further nextflow options]
# Do not add -resume: task inputs are dataset names only, so runs on the old truncated data would be reused.

profile=$1
shift

if [ "$profile" == "cibersortx" ]; then
    work_dir=/localscratch/omnideconv_work
elif [ "$profile" == "apptainer" ] || [ "$profile" == "autogenes" ]; then
    work_dir=/nfs/scratch/nf-core_work/omnideconv.benchmark
else
    echo "Usage: $0 <apptainer|autogenes|cibersortx> [further nextflow options]"
    exit 1
fi

nextflow -C /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/configs/nextflow_norm_morandini_ncbi.config run \
    /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline/main.nf \
    -profile $profile \
    -work-dir $work_dir \
    -with-trace /nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/reports_traces/morandini_ncbi-trace-$profile.csv "$@"
