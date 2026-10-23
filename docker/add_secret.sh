#!/bin/bash

# Script to create an image with the deconvolution method SECRET from the standard image.
# SECRET (https://github.com/anlingUA/SECRET) is not part of omnideconv. Only SECRET and the R
# package alabama, which it uses for the optimisation, are added; everything else stays as it is.
# Usage: ./add_secret.sh <source_directory> <new_name> <SECRET_commit>
#
# Example: ./add_secret.sh omnideconv_benchmark_2_0 omnideconv_benchmark_2.4 8ee3b53fcedcc33d04c79a5b5b1abca82824cf5e
#
# <source_directory> is an Apptainer sandbox directory (as created by docker_to_singularity.sh).
# Creates the directory <new_name> and the SIF file <new_name>.sif.
# Note: the SECRET repository does not state a licence.

set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <source_directory> <new_name> <SECRET_commit>"
    echo "Example: $0 omnideconv_benchmark_2_0 omnideconv_benchmark_2.4 8ee3b53fcedcc33d04c79a5b5b1abca82824cf5e"
    exit 1
fi

SOURCE_DIR="${1%/}"
OUTPUT_DIR="${2%/}"
COMMIT="$3"

if [ -e "$OUTPUT_DIR" ]; then
    echo "Error: $OUTPUT_DIR already exists"
    exit 1
fi

echo "Copying $SOURCE_DIR to $OUTPUT_DIR ..."
cp -a "$SOURCE_DIR" "$OUTPUT_DIR"

# a GitHub token that was used to build the source image must not be carried into new images
rm -f "$OUTPUT_DIR/root/.Renviron"

echo "Installing alabama and SECRET ($COMMIT) ..."
apptainer exec --writable --no-home "$OUTPUT_DIR" Rscript -e "
  # the image installs R packages through the system package manager (bspm), which needs root
  if (requireNamespace('bspm', quietly = TRUE)) bspm::disable()
  options(repos = c(CRAN = 'https://cloud.r-project.org'))
  install.packages('alabama')
  remotes::install_github('anlingUA/SECRET@$COMMIT', upgrade = 'never', dependencies = FALSE, build_vignettes = FALSE)
"

echo ""
echo "Installed in $OUTPUT_DIR:"
apptainer exec --no-home "$OUTPUT_DIR" Rscript -e "cat('alabama', as.character(packageVersion('alabama')), '\nSECRET', as.character(packageVersion('SECRET')), '\n')"

echo ""
echo "Building SIF file from directory..."
apptainer build "${OUTPUT_DIR}.sif" "$OUTPUT_DIR"

echo "Successfully created directory $OUTPUT_DIR and ${OUTPUT_DIR}.sif"
