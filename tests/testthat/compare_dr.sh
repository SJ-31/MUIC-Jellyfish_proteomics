#!/bin/bash

# models=("esm" "a2v" "prottrans")
models=("prottrans")
# tech_list=("umap" "tsne" "pca" "pcoa")
tech_list=("umap")
# metrics=("cosine")
metrics=("euclidean" "cosine")

for model in "${models[@]}"; do
    for tech in "${tech_list[@]}"; do
            for metric in "${metrics[@]}"; do
                Rscript ./dr_compare.r  \
                    -f "./output/protein_dr/${tech}_comparison_${model}"  \
                    -s C_indra \
                    --technique "$tech" \
                    --metric "$metric" \
                    --embedding_path "../nf-test-out/C_indra_${model}_embeddings/embeddings.hdf5" \
                    --dist_path "../nf-test-out/C_indra_${model}_embeddings/distances.hdf5" \
                    --r_source "$BIN"/R \
                    --python_source "$BIN" \
                    --combined_results ../../results/C_indra_A/1-First_pass/C_indra_all.tsv
            done
    done
done
