#!/bin/bash
# Run one workflow with several profiles, one after the other, in one launch directory with -resume.
# The profiles then share the tasks that do not depend on the method: the reference split and the
# simulation run once and are reused by the following profiles. The CIBERSORTx profile has its own
# work directory and repeats them; the seeds give it the same split and the same pseudo-bulks.
#
# Usage: run_profiles.sh <config> <workflow> [--check] [<profile> ...]
#
#   <config>     run config: a name from configs/ without '.config', or a path
#   <workflow>   as for run.sh
#   --check      only run the method bisque of the apptainer profile, to see that the run works.
#                Its split, pseudo-bulks and results are real and are reused by the full run.
#   <profile>    default: apptainer autogenes rectangle cibersortx
#
# Launch directory (Nextflow's cache and log): <benchmark>/launch/<config>

set -u

if [ $# -lt 2 ]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

scripts_dir=$(dirname "$(realpath "$0")")
pipeline_dir=$(dirname "$scripts_dir")
config=$1
workflow=$2
shift 2

if [ ! -f "$config" ]; then
    config=$pipeline_dir/configs/$config.config
fi
if [ ! -f "$config" ]; then
    echo "Config not found: $config"
    exit 1
fi
config=$(realpath "$config")

check=false
if [ $# -gt 0 ] && [ "$1" == "--check" ]; then
    check=true
    shift
fi
profiles=${*:-"apptainer autogenes rectangle cibersortx"}

launch_dir=$(dirname "$pipeline_dir")/launch/$(basename "$config" .config)
mkdir -p "$launch_dir"
cd "$launch_dir" || exit 1

if $check; then
    profiles=apptainer
    check_config=$launch_dir/$(basename "$config" .config)_check.config
    printf "// %s with one quick method\nincludeConfig '%s'\n\nprofiles {\n     apptainer {\n        params.method_list = [\"bisque\"]\n     }\n}\n" "$(basename "$config")" "$config" > "$check_config"
    config=$check_config
fi

status=0
for profile in $profiles; do
    echo "[$(date '+%F %T')] $(basename "$config" .config), profile $profile"
    "$scripts_dir/run.sh" "$config" "$profile" "$workflow" -resume > "nextflow_$profile.log" 2>&1
    code=$?
    echo "[$(date '+%F %T')] profile $profile finished with exit code $code (log: $launch_dir/nextflow_$profile.log)"
    failed=$(grep -c "Error executing process" "nextflow_$profile.log")
    if [ "$failed" -gt 0 ]; then
        echo "    $failed failed task(s):"
        grep "Error executing process" "nextflow_$profile.log" | sed "s/.*Error executing process > /      /" | sort -u
        status=1
    fi
    [ $code -eq 0 ] || status=1
done
exit $status
