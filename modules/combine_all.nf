process COMBINE_ALL {
    publishDir "$outdir", mode: "copy", pattern: "*.tsv"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(downloads)
    path(eggnog)
    path(interpro)
    tuple path(directlfq), path(directlfq_input), path(maxlfq)
    val(outdir)
    val(logdir)
    //

    output:
    path(output), emit: all
    path("combine_all.log")
    path("${params.pref}_taxonomy.tsv")
    //

    script:
    output = "${params.pref}_all.tsv"
    if (eggnog == "no_file") {
        eggnog = interpro = ""
    }
    if (file("${outdir}/${output}").exists()) {
     """
     cp ${outdir}/${output} .
     cp ${logdir}/combine_all.log .
     cp ${outdir}/${params.pref}_taxonomy.tsv .
     """
    } else {
    """
    Rscript ${params.bin}/R/combine_all.r \
        --eggnog $eggnog \
        --interpro $interpro \
        --downloads $downloads \
        --pfam2go $params.pfam2go \
        --interpro2go $params.interpro2go \
        --kegg2go $params.kegg2go \
        --ec2go $params.ec2go \
        --pfam_db $params.pfam_entries \
        --fdr $params.fdr \
        --denovo_org $params.species_spec \
        --pep_thresh $params.pep_thresh \
        --output "${params.pref}_all.tsv" \
        --go_path $params.go \
        --go_slim_path $params.go_slims \
        --directlfq_aqreformat $directlfq_input \
        --directlfq $directlfq \
        --maxlfq $maxlfq \
        --r_source ${params.bin}/R \
        --python_source ${params.bin}

    get_tax_data.py -i "${params.pref}_all.tsv" \
        -o  "${params.pref}_taxonomy.tsv" \
        -n  $params.ncbi_taxdump

    cat .command.log >> combine_all.log
    """
    }
    //
}
