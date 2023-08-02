#!/usr/bin/env sh

file_list=$1
outdir=$2
for format in {mzML,mzXML,mgf}
do 
    msconvert -f "$file_list" \
        -o "$outdir" \
        --"${format}" \
        --filter "peakPicking true MS2Deisotope true MS2Denoise true"
done

# Explore the options for peakPicking a bit more
