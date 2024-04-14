process EMBEDDINGS {
    publishDir "$outdir", mode: "copy", pattern: "*hdf5"
    publishDir "$params.logdir", mode: "copy", pattern: "*log"

    input:
    path(results_tsv)
    val(outdir)
    val(model)
    //

    output:
    path("embeddings.hdf5"), emit: embd
    path("distances.hdf5")
    path("*log")
    //

    script:
    def check = file("$outdir/embeddings.hdf5")
    if (check.exists()) {
        """
        cp "$outdir/embeddings.hdf5" .
        cp "$outdir/distances.hdf5" .
        cp "$params.logdir/${model}.log" .
        """
    } else {
        if (model == "prottrans") {
            """
            get_fasta.py -i ${results_tsv} -s seq -d ProteinId -o seqs.fasta

            source activate $params.protlm
            prott5_embedder.py --input seqs.fasta \
            --model "$params.prottrans" \
            --output embeddings.hdf5 \
            --per_protein 1

            get_distances.py -i embeddings.hdf5 \
                -o distances.hdf5
            conda deactivate
            cp .command.log "${model}".log
            """
        } else if (model == "esm") {
            """
            get_fasta.py -i ${results_tsv} -s seq -d ProteinId -o seqs.fasta

            source activate $params.esmfold
            mkdir output
            esm-extract "$params.esm" \
                seqs.fasta \
                output \
                --include mean

            get_distances.py -i output \
                -o distances.hdf5 \
                -w embeddings.hdf5
            conda deactivate
            cp .command.log "${model}".log
            """
        } else if (model == "a2v") {
            """
            Rscript $params.bin/R/analysis/prepare_embeddings.r \
                --r_source "$params.bin/R" \
                --python_source "$params.bin" \
                --sample_tsv $results_tsv \
                --embd_output embeddings.hdf5 \
                --dist_output distances.hdf5 \
                --embedding_path $params.go_embeddings

            cp .command.log "${model}".log
            """
        } else if (model == "a2v_go") {
            """
            Rscript $params.bin/R/analysis/prepare_embeddings.r \
                --r_source "$params.bin/R" \
                --python_source "$params.bin" \
                --sample_tsv $results_tsv \
                --embd_output embeddings.hdf5 \
                --dist_output distances.hdf5 \
                --embedding_path $params.go_embeddings \
                --go_only

            cp .command.log "${model}".log
            """
        }
    }
    //
}
