#!/usr/bin/bash

script=$(realpath ../../bin/R/ontologizer.r)
r_source=$(realpath ../../bin/R)
slims=$(realpath ../../data/reference/goslim_generic.obo)
data=$(realpath ../../data/reference/.go_texts)
go=$(realpath ../../data/reference/go.obo)

cd ../../results/C_indra/Analysis/Ontologizer

Rscript "$script" \
    --r_source "$r_source" \
    --results_path . \
    --go_slim_path "$slims"  \
    --go_tm_dir "$data" \
    --go_path  "$go" \
    --mode word_cloud
