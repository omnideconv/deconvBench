#!/bin/bash
# Build .sif files from the Apptainer sandboxes the cluster uses and push them to the GitHub
# container registry, where the CI (.github/workflows/ci.yml) pulls them from:
#   ghcr.io/omnideconv/deconvbench-images          standard-2.0, autogenes-1.3, rectangle-2.5 (public)
#   ghcr.io/omnideconv/deconvbench-images-secret   2.4 (private: the SECRET repository states no licence)
#
# Once before:  apptainer registry login --username <github user> oras://ghcr.io
#               (password: a GitHub token with the scope write:packages)
# A new package is private. The repository omnideconv/deconvBench needs read access to a private
# package (package settings, "Manage Actions access"); a public package needs an organisation owner
# to allow public packages.
#
# Usage: ./push_ci_images.sh <build_dir> [<image> ...]    images: standard autogenes rectangle secret

set -e
docker_dir=$(dirname "$(realpath "$0")")
build_dir=$1; shift
images=${*:-"standard autogenes rectangle secret"}
declare -A sandbox=([standard]=omnideconv_benchmark_2_0 [autogenes]=omnideconv_benchmark_1_3
                    [rectangle]=omnideconv_benchmark_2.5 [secret]=omnideconv_benchmark_2.4)
declare -A target=([standard]=deconvbench-images:standard-2.0 [autogenes]=deconvbench-images:autogenes-1.3
                   [rectangle]=deconvbench-images:rectangle-2.5 [secret]=deconvbench-images-secret:2.4)
mkdir -p "$build_dir"
# apptainer unpacks the sandbox before packing it; /tmp of the login node is too small for that
export APPTAINER_TMPDIR=${APPTAINER_TMPDIR:-$build_dir/tmp} APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$build_dir/cache}
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
for image in $images; do
    sif=$build_dir/$image.sif
    [ -f "$sif" ] || apptainer build "$sif" "$docker_dir/${sandbox[$image]}"
    apptainer push "$sif" "oras://ghcr.io/omnideconv/${target[$image]}"
    echo "pushed ${target[$image]}"
done
