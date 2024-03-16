#!/bin/bash

# tech_list=("umap" "tsne")
tech_list=("umap")

for tech in "${tech_list[@]}"; do
    TECH=$(echo $tech | tr '[:lower:]' '[:upper:]')
    for e_mode in "mean" "sum"; do
        for metric in "euclidean" "cosine"; do
            Rscript ./dr_compare.r  \
                -f ./output/figures/${tech}_comparison_${e_mode}  \
                -s C_indra \
                -u ../../data/protein_databases/comparison_taxa \
                --technique "$TECH" \
                --ontologizer_path ../nf-test-out/ontologizer/ \
                --uniprot_embeddings ../../data/protein_databases/uniprot_embeddings_${e_mode}.tsv \
                --metric "$metric" \
                --protein_embedding_mode ${e_mode} \
                --embeddings_path ../../data/reference/go_embedded.npz \
                --r_source "$BIN"/R \
                --python_source "$BIN" \
                --combined_results ../../results/C_indra/1-First_pass/C_indra_all.tsv
        done
    done
done
