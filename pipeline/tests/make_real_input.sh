#!/bin/bash
# Create the input directory for the real-data check (tests/conf_real/): links to the datasets on
# the norm cluster, in the layout of the fixture. The pipeline versions under test write their
# simulated datasets next to these links, not into the data directories.
#
# Usage: make_real_input.sh <input_dir>

set -eu

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

data=/nfs/data/omnideconv_benchmarking_clean/data
mkdir -p "$1/single_cell" "$1/bulks"

for sc in singleCell/lambrechts singleCell/maynard singleCellXGeneCensus/HaoCleaned singleCellXGeneCensus/HaoCleaned-sampled; do
    ln -sfn "$data/$sc" "$1/single_cell/$(basename "$sc")"
done
for bulk in PBMC/hoek PBMC/finotello PBMC/morandini Tumor/vanderbilt_lung; do
    ln -sfn "$data/$bulk" "$1/bulks/$(basename "$bulk")"
done

ls -l "$1/single_cell" "$1/bulks"
