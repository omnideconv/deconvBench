#!/bin/bash
# Run every workflow of the pipeline on the regression fixture.
#
# Usage: run_fixture.sh <fixture_dir> <run_dir>
#
#   <fixture_dir>  output of make_fixture.R
#   <run_dir>      receives input/ (per-workflow copy of the fixture), results/, launch/ (logs, traces)
#                  has to be below /nfs/data/omnideconv_benchmarking_clean (mounted into the CIBERSORTx container)
#
# Environment:
#   CODE_DIR    pipeline directory whose main.nf is tested (default: the one this script is in)
#   PROFILES    default: "apptainer autogenes rectangle cibersortx"
#   WORKFLOWS   default: all eight
#   NF_ARGS     additional nextflow run options, e.g. "-resume"
#
# Profiles run one after the other, as in production, so that later profiles see the simulations
# of earlier ones. Within a Slurm profile the workflows run in parallel.

set -u

if [ $# -ne 2 ]; then
    echo "Usage: $0 <fixture_dir> <run_dir>"
    exit 1
fi

fixture_dir=$(realpath "$1")
mkdir -p "$2"
run_dir=$(realpath "$2")

tests_dir=$(dirname "$(realpath "$0")")
code_dir=$(realpath "${CODE_DIR:-$tests_dir/..}")
profiles=${PROFILES:-"apptainer autogenes rectangle cibersortx"}
workflows=${WORKFLOWS:-"main simulation subsampling impact_missing_cell_types simulation_spillover simulation_unknown_content simulation_resolution_analysis simulation_impact_technology"}
nf_args=${NF_ARGS:-}

run_workflow() {
    local profile=$1
    local wf=$2
    local launch_dir=$run_dir/launch/$profile-$wf
    local work_dir=/nfs/scratch/nf-core_work/omnideconv.refactor_test/$(basename "$run_dir")/$profile-$wf
    if [ "$profile" == "cibersortx" ]; then
        work_dir=/localscratch/omnideconv_work/refactor_test/$(basename "$run_dir")/$wf
    fi
    local entry=""
    if [ "$wf" != "main" ]; then
        entry="-entry $wf"
    fi

    mkdir -p "$launch_dir" "$run_dir/results/$wf"
    cd "$launch_dir" || return 1
    TEST_INPUT=$run_dir/input/$wf TEST_OUT=$run_dir/results/$wf NXF_ANSI_LOG=false \
        nextflow -C "$tests_dir/conf/$wf.config" run "$code_dir/main.nf" \
        -profile "$profile" $entry \
        -work-dir "$work_dir" \
        -with-trace "$launch_dir/trace.csv" $nf_args > nextflow.out 2>&1
    echo $? > exit_status
}

for wf in $workflows; do
    if [ ! -d "$run_dir/input/$wf" ]; then
        mkdir -p "$run_dir/input"
        cp -r "$fixture_dir" "$run_dir/input/$wf"
        mkdir -p "$run_dir/input/$wf/preprocess"
    fi
done

for profile in $profiles; do
    echo "[$(date +%H:%M:%S)] profile $profile"
    for wf in $workflows; do
        if [ "$profile" == "cibersortx" ]; then
            run_workflow "$profile" "$wf"
        else
            run_workflow "$profile" "$wf" &
        fi
    done
    wait
done

echo
printf '%-12s %-34s %s\n' profile workflow exit
status=0
for profile in $profiles; do
    for wf in $workflows; do
        code=$(cat "$run_dir/launch/$profile-$wf/exit_status" 2>/dev/null || echo missing)
        printf '%-12s %-34s %s\n' "$profile" "$wf" "$code"
        [ "$code" == "0" ] || status=1
    done
done
exit $status
