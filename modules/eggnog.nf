process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"
    memory "35 GB"
    errorStrategy 'ignore'

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

    shell:
    check = file("${outdir}/${params.pref}.emapper.seed_orthologs")
    if (check.exists()) {
        '''
        cp !{outdir}/!{params.pref}.emapper* .
        '''
    } else {
        '''
        clean_fasta.py -i !unknown_fasta -o query.fasta
        if [[ -e !{params.saved}/eggnog_annotations.tsv ]]; then
            helpers.py -t get_saved \
                --save_type eggnog \
                -i query.fasta \
                --saved !{params.saved} \
                --seq_header_mapping !{results}/Databases/seq-header_mappings.tsv
            rm query.fasta; mv new_query.fasta query.fasta
        fi

        export EGGNOG_DATA_DIR="!params.eggnog_data_dir"
        emapper.py -i query.fasta \
            --output !{params.pref} \
            -m mmseqs \
            --report_orthologs \
            --no_file_comments

        if [[ -e saved_annotations.tsv ]]; then
            for i in "annotations" "orthologs" "seed_orthologs"; do
                cat saved_${i}.tsv >> !{params.pref}.emapper.${i} saved_${i}.tsv
            done
        fi
        '''
    }
}

process SORT_EGGNOG {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(eggnog_annotations), path(eggnog_seeds),\
        path(unknown_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}_eggnog.fasta"), emit: fasta
    path("${params.pref}_eggnog_unmatched.tsv"), emit: unmatched
    path("${params.pref}_eggnog_matched.tsv"), emit: matched
    //

    script:
    def check = file("${outdir}/${params.pref}_eggnog_matched.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/*${params.pref}_eggnog* .
        """
    } else {
        """
        Rscript $params.bin/R/sort_eggnog.r \
            -u ${params.pref}_eggnog_unmatched.tsv \
            -f ${params.pref}_eggnog.fasta \
            --blast $unknown_tsv \
            --annotations $eggnog_annotations \
            --seeds $eggnog_seeds \
            --output ${params.pref}_eggnog_matched.tsv
        """
    }
    //
}
