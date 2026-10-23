#!/bin/bash
# Build .sif files from the Apptainer sandboxes the cluster uses and push them to the GitHub
# container registry, where the CI (.github/workflows/ci.yml) pulls them from.
#
# Once before:  apptainer registry login --username <github user> oras://ghcr.io
#               (password: a GitHub token with the scope write:packages)
# Afterwards:   make the package ghcr.io/omnideconv/deconvbench-images public on GitHub.
#
# Usage: ./push_ci_images.sh <build_dir> [<tag> ...]    tags: standard-2.0 autogenes-1.3 rectangle-2.5 secret-2.4

set -e
docker_dir=$(dirname "$(realpath "$0")")
build_dir=$1; shift
tags=${*:-"standard-2.0 autogenes-1.3 rectangle-2.5 secret-2.4"}
declare -A sandbox=([standard-2.0]=omnideconv_benchmark_2_0 [autogenes-1.3]=omnideconv_benchmark_1_3
                    [rectangle-2.5]=omnideconv_benchmark_2.5 [secret-2.4]=omnideconv_benchmark_2.4)
mkdir -p "$build_dir"
# apptainer unpacks the sandbox before packing it; /tmp of the login node is too small for that
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-$build_dir/tmp} APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$build_dir/cache}
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
for tag in $tags; do
    sif=$build_dir/$tag.sif
    [ -f "$sif" ] || apptainer build "$sif" "$docker_dir/${sandbox[$tag]}"
    apptainer push "$sif" "oras://ghcr.io/omnideconv/deconvbench-images:$tag"
    echo "pushed $tag"
done
