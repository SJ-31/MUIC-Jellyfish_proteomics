#!/usr/bin/bash
protlm="/home/shannc/anaconda3/envs/protlm"
data="/mnt/data/shannc/nf/data"
model="$data/models/prot_t5_xl_uniref50"
bin="/mnt/data/shannc/nf/bin"

results="$1"
seq_col="$2"
id_col="$3"
"$bin/get_fasta.py" -i "$results" -s "$seq_col" -d "$id_col" -o seqs.fasta

source activate "$protlm"
prott5_embedder.py --input seqs.fasta \
    --model "$model" \
    --output embeddings.hdf5 \
    --per_protein 1

"$bin/get_distances.py" -i embeddings.hdf5 \
    -o distances.hdf5
conda deactivate
