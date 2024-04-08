#!/usr/bin/env bash

pt=/home/shannc/workflow/data/models/prot_t5_xl_uniref50
esm=/home/shannc/workflow/data/models/esm2_t36_3B_UR50D.pt
uniprot_dir="/home/shannc/workflow/data/protein_databases/comparison_taxa"
outdir="/home/shannc/workflow/data/protein_databases/uniprot_embeddings"

cd "$outdir" || exit
for file in $(find "$uniprot_dir" -name "*reviewed.tsv"); do
    out=$(echo "${file}_seqs" | sed 's;.*/\(.*\).tsv;\1;')
    out="${out}.temp"
    "$BIN"/get_fasta.py -i "$file" -s Sequence -d Entry -o "$out"
    echo "line count: $(wc -l "$file"), num seqs = $(seqkit stat "$out" -T | cut -f 4 | tail -n 1)"
done

cat ./*temp | seqkit seq -m 10 > uniprot_sequences.fasta
cat uniprot_sequences.fasta | grep ">" | grep -v "^>"

source activate protlm
prott5_embedder.py --input uniprot_sequences.fasta \
    --model "$pt" \
    --output pt_embeddings.hdf5 \
    --per_protein 1
conda deactivate

source activate esmfold
mkdir "${outdir}"/esm_out
esm-extract "$esm" \
    uniprot_sequences.fasta \
    esm_out \
    --include mean
conda deactivate
