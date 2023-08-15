#!/usr/bin/env sh

for t in tide*percolator*
    do
    rename="$(echo "$t" | sed -e 's/\./_/g' -e 's/_txt/\.tsv/' -e 's/_target//')"
    mv $t $rename
done
