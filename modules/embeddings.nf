process EMBEDDINGS {
    publishDir "$outdir", mode: ""

    input:
    path(results_tsv)
    val(outdir)
    val(model)
    //

    output:
    path("embeddings.hdf5")
    path("distances.hdf5")
    //

    script:
    def check = file("$outdir/embeddings.hdf5")
    if (check.exists()) {
        """
        cp "$outdir/embeddings.hdf5" .
        cp "$outdir/distances.hdf5" .
        """
    } else {
        if (model == "prottrans") {
            """
            get_fasta.py -i ${results_tsv} -s seq -d ProteinId -o seqs.fasta

            source activate /home/shannc/anaconda3/envs/protlm
            prott5_embedder.py --input seqs.fasta \
            --model "$params.prottrans" \
            --output embeddings.hdf5 \
            --per_protein 1

            get_distances.py -i embeddings.hdf5 \
                -o distances.hdf5
            conda deactivate
            """
        } else if (model == "esm") {
            """
            get_fasta.py -i ${results_tsv} -s seq -d ProteinId -o seqs.fasta

            source activate /home/shannc/anaconda3/envs/esmfold
            mkdir output
            esm-extract "$params.esm" \
                seqs.fasta \
                output \
                --include mean

            get_distances.py -i output \
                -o distances.hdf5 \
                -w embeddings.hdf5
            conda deactivate
            """
        }
    }
    //
}
