process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"

    input:
    tuple path(unknown_fasta), path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}.emapper.seed_orthologs"), emit: seed
    path("${params.pref}.emapper.orthologs"), emit: orthologs
    path("${params.pref}.emapper.annotations"), emit: annotations
    tuple path("${params.pref}.emapper.annotations"), \
        path("${params.pref}.emapper.seed_orthologs"), \
        path(unknown_tsv), emit: unmatched
    //

    script:
    check = file("${outdir}/${params.pref}.emapper.seed_orthologs")
    if (check.exists()) {
        """
        cp ${outdir}/${params.pref}.emapper* .
        """
    } else {
        """
        export EGGNOG_DATA_DIR="$params.eggnog_data_dir"
        emapper.py -i $unknown_fasta \
            --output ${params.pref} \
            -m mmseqs \
            --report_orthologs
        """
    }
}

process SORT_EGGNOG {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(eggnog_annotations), path(eggnog_seeds),\
        path(unknown_tsv)
    path(unmatched_peptides)
    val(outdir)
    //

    output:
    path("${params.pref}_eggnog.fasta"), emit: fasta
    path("${params.pref}_eggnog_unmatched.tsv"), emit: unmatched
    path("${params.pref}_eggnog_matched.tsv"), emit: matched
    //

    script:
    def check = file("${params.pref}_eggnog_matched.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/*${params.pref}_eggnog* .
        cp ${outdir}/eggnog_*${params.pref}* .
        """
    } else {
        """
        Rscript $params.bin/sort_eggnog.r \
            -u ${params.pref}_eggnog_unmatched.tsv \
            -f ${params.pref}_eggnog.fasta \
            --blast $unknown_tsv \
            --annotations $eggnog_annotations \
            --seeds $eggnog_seeds \
            --output temp.tsv

        eggnog_seq.py -e $params.eggnog_db_fasta \
            -i temp.tsv \
            -o ${params.pref}_eggnog_matched.tsv \
        """
    }
    //
}
