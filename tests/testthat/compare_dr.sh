#!/bin/bash

models=("esm" "a2v" "prottrans")
tech_list=("umap" "tsne")
tech_list2=("pca" "pcoa")

for model in "${models[@]}"; do
    for tech in "${tech_list[@]}"; do
        TECH=$(echo $tech | tr '[:lower:]' '[:upper:]')
            for metric in "euclidean" "cosine"; do
                Rscript ./dr_compare.r  \
                    -f "./output/protein_dr/${tech}_comparison"  \
                    -s C_indra \
                    --technique "$TECH" \
                    --metric "$metric" \
                    --embedding_path "../nf-test-out/C_indra_${model}_embeddings/embeddings.hdf5" \
                    --dist_path "../nf-test-out/C_indra_${model}_embeddings/distances.hdf5" \
                    --r_source "$BIN"/R \
                    --python_source "$BIN" \
                    --combined_results ../../results/C_indra/1-First_pass/C_indra_all.tsv
            done
    done
    # for tech in "${tech_list2[@]}"; do
    #     # Rscript "$BIN"/R/


    # done
done
