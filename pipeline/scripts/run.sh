#!/bin/bash
# Start a workflow of the benchmark on the DaisyBio cluster.
#
# Usage: run.sh <config> <profile> [<workflow>] [further nextflow options]
#
#   <config>     experiment config: a name from configs/ without '.config', or a path
#   <profile>    apptainer | autogenes | rectangle | cibersortx  (see conf/norm.config)
#   <workflow>   main (default) | deconv_only | simulation | subsampling | impact_missing_cell_types |
#                simulation_spillover | simulation_unknown_content | simulation_resolution_analysis |
#                simulation_impact_technology
#
# Examples:
#   run.sh nextflow_norm_morandini_ncbi apptainer
#   run.sh nextflow_norm_resolution rectangle simulation_resolution_analysis -resume
#
# A trace file is written to reports_traces/<config>-<workflow>-<profile>-trace.csv.
# -resume is safe to use: tasks are re-run when their input files change.

set -u

if [ $# -lt 2 ]; then
    sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

pipeline_dir=$(dirname "$(dirname "$(realpath "$0")")")
config=$1
profile=$2
shift 2

if [ ! -f "$config" ]; then
    config=$pipeline_dir/configs/$config.config
fi
if [ ! -f "$config" ]; then
    echo "Config not found: $config"
    exit 1
fi

workflow=main
if [ $# -gt 0 ] && [[ "$1" != -* ]]; then
    workflow=$1
    shift
fi
entry=()
if [ "$workflow" != "main" ]; then
    entry=(-entry "$workflow")
fi

# CIBERSORTx runs in Docker on this node and needs a work directory on a local disk
if [ "$profile" == "cibersortx" ]; then
    work_dir=/localscratch/omnideconv_work
else
    work_dir=/nfs/scratch/nf-core_work/omnideconv.benchmark
fi

# With Java 21 on this cluster, copying a file larger than about 2 MB to an NFS share never returns,
# which blocks Nextflow when it publishes results. Java 17 is not affected.
export NXF_JAVA_HOME=${NXF_JAVA_HOME:-/usr/lib/jvm/java-17-openjdk-amd64}

trace_file=$(dirname "$pipeline_dir")/reports_traces/$(basename "$config" .config)-$workflow-$profile-trace.csv

nextflow -C "$config" run "$pipeline_dir/main.nf" \
    -profile "$profile" "${entry[@]}" \
    -work-dir "$work_dir" \
    -with-trace "$trace_file" "$@"
