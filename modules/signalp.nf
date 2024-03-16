process SIGNALP {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/signalp"

    input:
    path(unmatched_fasta)
    val(outdir)
    //

    output:
    path("processed_entries.fasta")
    path("prediction_results.txt")
    path("*.gff3")

    //

    script:
    def check = file("${outdir}/processed_entries.fasta")
    if (check.exists()) {
        """
        mv ${outdir}/processed_entries.fasta .
        mv ${outdir}/prediction_results.txt .
        mv ${outdir}/*.gff3 .
        """
    } else {
        """
        mkdir results
        seqkit seq -m 10 -g $unmatched_fasta > filtered.fasta
        signalp6 -ff filtered.fasta \
            -org euk \
            --output_dir results \
            --model_dir ${params.signalp6model}
        cp results/prediction_results.txt .
        cp results/processed_entries.fasta .
        cp results/*.gff3 .
        """
    }

    //
}
