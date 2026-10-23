#!/bin/bash
# TCGA-BRCA deconvolved with wu (configs/tcga_brca.config): eight omnideconv methods, rectangle, secret.
# Usage: run_tcga_brca.sh [--check] [<profile> ...]
scripts_dir=$(dirname "$(realpath "$0")")
check=()
if [ "${1:-}" == "--check" ]; then check=(--check); shift; fi
"$scripts_dir/run_profiles.sh" tcga_brca main "${check[@]}" ${*:-apptainer autogenes rectangle secret cibersortx}
