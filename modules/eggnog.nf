process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"

    input:
    tuple path(unknown_fasta), path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${name}.emapper.seed_orthologs"), emit: seed
    path("${name}.emapper.orthologs"), emit: orthologs
    path("${name}.emapper.annotations"), emit: annotations
    tuple val("${name}"), path("${name}.emapper.annotations"), \
        path("${name}.emapper.seed_orthologs"), \
        path(unknown_tsv), emit: unmatched
    //

    script:
    name = unknown_fasta.baseName.replaceFirst(/_unmatched/, "")
    check = file("${outdir}/${name}.emapper.seed_orthologs")
    if (check.exists()) {
        """
        cp ${outdir}/${name}.emapper* .
        """
    } else {
        """
        export EGGNOG_DATA_DIR="$params.eggnog_data_dir"
        emapper.py -i $unknown_fasta \
            --output ${name} \
            --report_orthologs
        """
    }
}

process SORT_EGGNOG {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(run_name), path(eggnog_annotations), path(eggnog_seeds),\
        path(unknown_tsv)
    path(unmatched_peptides)
    val(outdir)
    //

    output:
    path("${run_name}_eggnog.fasta"), emit: fasta
    path("${run_name}_eggnog_unmatched.tsv"), emit: unmatched
    tuple val(run_name), path("eggnog_anno-${run_name}.tsv"), path("eggnog_meta-${run_name}.tsv"), emit: matched
    //

    script:
    def check = file("${outdir}/eggnog_anno-${run_name}.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/*${run_name}_eggnog* .
        cp ${outdir}/eggnog_*${run_name}* .
        """
    } else {
        """
        Rscript $params.bin/sort_eggnog.r \
            -u ${run_name}_eggnog_unmatched.tsv \
            -f ${run_name}_eggnog.fasta \
            --blast $unknown_tsv \
            --annotations $eggnog_annotations \
            --seeds $eggnog_seeds \
            --output_meta meta_temp \
            --output_anno "eggnog_anno-${run_name}.tsv"

        eggnog_seq.py -e $params.eggnog_db_fasta \
            -m meta_temp \
            -a eggnog_anno-${run_name}.tsv \
            -o eggnog_meta-${run_name}.tsv
        """
    }
    //
}
