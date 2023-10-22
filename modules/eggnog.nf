process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"

    input:
    tuple path(unknown_fasta), path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${unknown_fasta.baseName}.emapper.seed_orthologs")
    path("${unknown_fasta.baseName}.emapper.orthologs")
    path("${unknown_fasta.baseName}.emapper.annotations")
    //

    script:
    """
    export EGGNOG_DATA_DIR="$params.eggnog_data_dir"
    emapper.py -i $unknown_fasta \
        --output ${unknown_fasta.baseName} \
        --report_orthologs

    unmatched_eggnog.r -f ${unknown_fasta.baseName}_eggnog.fasta \
        -b $unknown_tsv \
        -a ${unknown_fasta.baseName}.emapper.annotations
    """
    //
}
