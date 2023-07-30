#!/usr/bin/env sh

file_list=$1
outdir=$2
msconvert -f "$file_list" \
    -o "$outdir" \
    --mzML \
    --mzXML \
    --filter "peakPicking vendor MS2Deisotope"

# Explore the options for peakPicking a bit more
