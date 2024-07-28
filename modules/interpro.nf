process INTERPROSCAN {
    publishDir "$outdir", mode: "copy"
    memory "35 GB"

    input:
    path(unknown_fasta)
    val(outdir)
    //

    output:
    path("${params.pref}_interpro.tsv")
    //

    shell:
    check = file("${outdir}/${params.pref}_interpro.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/!{params.pref}_interpro.tsv .
        '''
    } else {
        '''
        clean_fasta.py -i !{unknown_fasta} -o query.fasta
        header="query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways"
        if [[ -e !{params.saved}/interpro.tsv ]]; then
            helpers.py -t get_saved \
                --save_type interpro \
                -i query.fasta \
                --saved !{params.saved}/interpro.tsv \
                --seq_header_mapping !{results}/Databases/seq-header_mappings.tsv
            rm query.fasta; mv new_query.fasta query.fasta

        interproscan.sh -i query.fasta \
            -f tsv \
            -goterms \
            -pa \
            -b temp

        if [[ -e saved_interpro.tsv ]]; then
            cat saved_interpro.tsv >> temp.tsv
        fi

        echo -e "$header" > header.txt
        cat header.txt temp.tsv > !{params.pref}_interpro.tsv
        '''
    }
    //
}

process SORT_INTERPRO {
    publishDir "$outdir", mode: "copy"
    publishDir "$outdir_unmatched", mode: "copy", pattern: "remaining_unmatched*"

    input:
    path(interpro_results)
    path(eggnog_unmatched)
    val(outdir)
    val(outdir_unmatched)
    //

    output:
    path("${params.pref}_interpro_matched.tsv"), emit: matched
    path("remaining_unmatched-${params.pref}.tsv"), emit: unmatched
    path("remaining_unmatched-${params.pref}.fasta"), emit: unmatched_fasta
    //

    script:
    """
    Rscript $params.bin/R/sort_interpro.r \
        -i $interpro_results \
        -u $eggnog_unmatched \
        -o ${params.pref}_interpro_matched.tsv \
        -f remaining_unmatched-${params.pref}.tsv \
        -a remaining_unmatched-${params.pref}.fasta
    """
    //
}
