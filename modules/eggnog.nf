process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"

    input:
    tuple path(unknown_fasta), path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${unknown_fasta.baseName}.emapper.seed_orthologs"), emit: seed
    path("${unknown_fasta.baseName}.emapper.orthologs"), emit: orthologs
    path("${unknown_fasta.baseName}.emapper.annotations"), emit: annotations
    tuple val("${unknown_fasta.baseName}"), path("${unknown_fasta.baseName}.emapper.annotations"), path(unknown_tsv), emit: unmatched
    //

    script:
    """
    export EGGNOG_DATA_DIR="$params.eggnog_data_dir"
    emapper.py -i $unknown_fasta \
        --output ${unknown_fasta.baseName} \
        --report_orthologs
    """
    //
}

process SORT_EGGNOG {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(run_name), path(eggnog_annotations), path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${run_name}_eggnog.fasta")
    //

    script:
    """
    Rscript $params.bin/unmatched_eggnog.r -f ${run_name}_eggnog.fasta \
        -b $unknown_tsv \
        -a $eggnog_annotations
    """
    //
}
