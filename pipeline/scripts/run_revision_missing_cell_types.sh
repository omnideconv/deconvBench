#!/bin/bash
# Revision: missing-cell-types analysis with HaoCleaned-sampled as reference (no split).
# Eight omnideconv methods and rectangle (configs/revision_missing_cell_types.config).
# Usage: run_revision_missing_cell_types.sh [--check] [<profile> ...]
scripts_dir=$(dirname "$(realpath "$0")")
check=()
if [ "${1:-}" == "--check" ]; then check=(--check); shift; fi
"$scripts_dir/run_profiles.sh" revision_missing_cell_types impact_missing_cell_types "${check[@]}" ${*:-apptainer autogenes rectangle cibersortx}
