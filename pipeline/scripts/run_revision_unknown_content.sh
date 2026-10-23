#!/bin/bash
# Revision: unknown-content analysis on lambrechts with the reference split by donor.
# Eight omnideconv methods, rectangle, secret and music-hc (configs/revision_unknown_content.config).
# Usage: run_revision_unknown_content.sh [--check] [<profile> ...]
scripts_dir=$(dirname "$(realpath "$0")")
check=()
if [ "${1:-}" == "--check" ]; then check=(--check); shift; fi
"$scripts_dir/run_profiles.sh" revision_unknown_content simulation_unknown_content "${check[@]}" ${*:-apptainer autogenes rectangle secret cibersortx}
