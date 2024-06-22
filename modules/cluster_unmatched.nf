process CLUSTER_UNMATCHED {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"
    errorStrategy 'ignore'

    input:
    path(remaining_unmatched)
    val(outdir)
    //

    output:
    path("unmatched_clustered.fasta"), emit: fasta
    path("unmatched_clustered.tsv")
    //

    shell:
    template 'cluster_unmatched.sh'
    //
}
