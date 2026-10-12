#!/bin/bash
# Revision: resolution analysis on lambrechts and wu with the reference split by donor.
# Eight omnideconv methods, rectangle and music-hc (configs/revision_resolution_<dataset>.config).
# Usage: run_revision_resolution.sh [--check] [<profile> ...]
scripts_dir=$(dirname "$(realpath "$0")")
check=()
if [ "${1:-}" == "--check" ]; then check=(--check); shift; fi
status=0
for dataset in lambrechts wu; do
    "$scripts_dir/run_profiles.sh" revision_resolution_$dataset simulation_resolution_analysis "${check[@]}" ${*:-apptainer autogenes rectangle cibersortx} || status=1
done
exit $status
