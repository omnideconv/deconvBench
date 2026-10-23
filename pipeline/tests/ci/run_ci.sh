#!/bin/bash
# Run one CI suite and check its results.
#
# Usage: run_ci.sh <suite> <run_dir>
#   workflows        all eight workflows with MuSiC, without split, split by cell and (two workflows)
#                    split by donor, then check_results.R and unit_tests.R
#   standard         main workflow with all methods of the standard image
#   autogenes | rectangle | secret   main workflow with the method of that image
#
# Environment: NEXTFLOW (launcher, default: nextflow), IMG_STANDARD, IMG_AUTOGENES, IMG_RECTANGLE, IMG_SECRET (images), CI_BIND (directories
# to bind into the containers; default: the repository and <run_dir>), NXF_ARGS (extra options).

set -u
suite=$1
mkdir -p "$2"
run_dir=$(realpath "$2")
ci_dir=$(dirname "$(realpath "$0")")
pipeline_dir=$(realpath "$ci_dir/../..")
repo_dir=$(dirname "$pipeline_dir")
export CI_BIND=${CI_BIND:-"$repo_dir,$run_dir"}

# <id> <workflow> [nextflow params]
run() {
    local id=$1 wf=$2; shift 2
    local entry=()
    [ "$wf" != "main" ] && entry=(-entry "$wf")
    # impact_missing_cell_types needs a ground truth for every bulk dataset
    [ "$wf" == "impact_missing_cell_types" ] && entry+=(-params-file "$ci_dir/params_missing_cell_types.json")
    mkdir -p "$run_dir/$id"
    rm -rf "$run_dir/$id/input"; cp -r "$repo_dir/data/ci" "$run_dir/$id/input"
    echo "[$(date +%T)] $id"
    (cd "$run_dir/$id" && CI_INPUT=$run_dir/$id/input CI_OUT=$run_dir/$id/results NXF_ANSI_LOG=false \
        "${NEXTFLOW:-nextflow}" -q -C "$ci_dir/ci.config" run "$pipeline_dir/main.nf" "${entry[@]}" ${PROFILE:+-profile $PROFILE} \
        -work-dir "$run_dir/work" ${NXF_ARGS:-} "$@" > nextflow.log 2>&1)
    local code=$?
    if [ $code -ne 0 ]; then echo "FAILED: $id (exit $code)"; tail -30 "$run_dir/$id/nextflow.log"; return 1; fi
}

status=0
case $suite in
  workflows)
    for wf in main simulation subsampling impact_missing_cell_types simulation_spillover simulation_unknown_content simulation_resolution_analysis simulation_impact_technology; do
        run "$wf" "$wf" || status=1
    done
    for wf in simulation subsampling impact_missing_cell_types simulation_spillover simulation_unknown_content simulation_resolution_analysis simulation_impact_technology; do
        run "$wf-cellsplit" "$wf" --split_reference true --split_by cell || status=1
    done
    for wf in simulation_unknown_content simulation_resolution_analysis; do
        run "$wf-donorsplit" "$wf" --split_reference true --split_by donor || status=1
    done
    ;;
  standard)  PROFILE=standard_methods run main-standard main || status=1 ;;
  autogenes|rectangle|secret)  PROFILE=$suite run "main-$suite" main || status=1 ;;
  *) echo "unknown suite: $suite"; exit 2 ;;
esac

# checks (base R only) run in the image of the suite, so that a method job needs a single image
case $suite in
  autogenes) check_image=$IMG_AUTOGENES ;;
  rectangle) check_image=$IMG_RECTANGLE ;;
  secret)    check_image=$IMG_SECRET ;;
  *)         check_image=$IMG_STANDARD ;;
esac
echo "[$(date +%T)] checks"
apptainer exec --no-home --cleanenv --bind "$CI_BIND" "$check_image" \
    Rscript "$ci_dir/check_results.R" "$suite" "$run_dir" "$pipeline_dir" || status=1
# unit tests: in the standard image with the workflows, in the SECRET image for SECRET's reference
case $suite in
  workflows) unit_image=$IMG_STANDARD ;;
  secret)    unit_image=$IMG_SECRET ;;
  *)         unit_image= ;;
esac
if [ -n "$unit_image" ]; then
    apptainer exec --no-home --cleanenv --bind "$CI_BIND" "$unit_image" \
        Rscript "$ci_dir/unit_tests.R" "$repo_dir" "$run_dir/unit" || status=1
fi
exit $status
