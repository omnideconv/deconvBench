#!/bin/bash

# Script to convert Docker image to Singularity SIF file and extract to folder
# Usage: ./docker_to_singularity.sh <docker_image_name> [output_name]
#
# Example: ./docker_to_singularity.sh omnideconv_benchmark:2.2
#          ./docker_to_singularity.sh omnideconv_benchmark:2.2 omnideconv_benchmark_2_2

set -e

# Check if docker image name is provided
if [ -z "$1" ]; then
    echo "Error: Docker image name is required"
    echo "Usage: $0 <docker_image_name> [output_name]"
    echo "Example: $0 omnideconv_benchmark:2.2"
    echo "Example: $0 omnideconv_benchmark:2.2 omnideconv_benchmark_2_2"
    exit 1
fi

DOCKER_IMAGE="$1"

# Determine output name
if [ -z "$2" ]; then
    # Convert docker image name to valid filename (replace : and / with _)
    OUTPUT_NAME=$(echo "$DOCKER_IMAGE" | sed 's/[:\/]/_/g')
else
    OUTPUT_NAME="$2"
fi

SIF_FILE="${OUTPUT_NAME}.sif"
OUTPUT_DIR="${OUTPUT_NAME}"

echo "Converting Docker image: $DOCKER_IMAGE"
echo "Output SIF file: $SIF_FILE"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Check if apptainer is available
if ! command -v apptainer &> /dev/null; then
    echo "Error: apptainer command not found. Please install Apptainer."
    exit 1
fi

# Build SIF file from Docker image
echo "Building SIF file from Docker image..."
apptainer build "$SIF_FILE" "docker-daemon://$DOCKER_IMAGE"

if [ $? -ne 0 ]; then
    echo "Error: Failed to build SIF file"
    exit 1
fi

echo "Successfully created $SIF_FILE"
echo ""

# Build sandbox (directory) from SIF file
echo "Extracting SIF file to directory..."
apptainer build --sandbox "$OUTPUT_DIR" "$SIF_FILE"

if [ $? -ne 0 ]; then
    echo "Error: Failed to extract SIF file to directory"
    exit 1
fi

echo "Successfully created directory: $OUTPUT_DIR"
echo ""
echo "Conversion complete!"
echo "  SIF file: $SIF_FILE"
echo "  Directory: $OUTPUT_DIR"
