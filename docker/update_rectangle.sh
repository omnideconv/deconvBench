#!/bin/bash

# Script to create a new rectangle image from an existing one with another version of rectanglepy.
# Only rectanglepy is changed; R, omnideconv, SimBu and all other Python packages stay as they are.
# Usage: ./update_rectangle.sh <source_directory> <new_name> <rectanglepy_version>
#
# Example: ./update_rectangle.sh omnideconv_benchmark_2.2 omnideconv_benchmark_2.3 1.6.0
#
# <source_directory> is an Apptainer sandbox directory (as created by docker_to_singularity.sh).
# Creates the directory <new_name> and the SIF file <new_name>.sif.
# To build an image from scratch instead, set the version in docker_rectangle/environment.yml,
# build the Docker image and convert it with docker_to_singularity.sh.

set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <source_directory> <new_name> <rectanglepy_version>"
    echo "Example: $0 omnideconv_benchmark_2.2 omnideconv_benchmark_2.3 1.6.0"
    exit 1
fi

SOURCE_DIR="${1%/}"
OUTPUT_DIR="${2%/}"
VERSION="$3"
PYTHON=/root/.local/share/r-miniconda/envs/r-omnideconv/bin/python

if [ -e "$OUTPUT_DIR" ]; then
    echo "Error: $OUTPUT_DIR already exists"
    exit 1
fi

echo "Copying $SOURCE_DIR to $OUTPUT_DIR ..."
cp -a "$SOURCE_DIR" "$OUTPUT_DIR"

# the GitHub token that was used to build the source image must not be carried into new images
rm -f "$OUTPUT_DIR/root/.Renviron"

echo "Installing rectanglepy $VERSION ..."
apptainer exec --writable --no-home "$OUTPUT_DIR" "$PYTHON" -m pip install --no-cache-dir "rectanglepy==$VERSION"

echo ""
echo "Python packages of rectangle in $OUTPUT_DIR:"
apptainer exec --no-home "$OUTPUT_DIR" "$PYTHON" -m pip list 2>/dev/null | grep -iE "^(rectanglepy|anndata|numpy|scipy|pandas|pydeseq2|osqp|statsmodels) "

echo ""
echo "Building SIF file from directory..."
apptainer build "${OUTPUT_DIR}.sif" "$OUTPUT_DIR"

echo "Successfully created directory $OUTPUT_DIR and ${OUTPUT_DIR}.sif"
