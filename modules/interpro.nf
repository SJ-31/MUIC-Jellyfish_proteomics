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
        clean_fasta.py -i !{unknown.fasta} -o query.fasta
        header="query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways"
        interproscan.sh -i query.fasta \
            -f tsv \
            -goterms \
            -pa \
            -b temp
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
    path("remaining_unmatched-${params.pref}.tsv")
    //

    script:
    """
    Rscript $params.bin/R/sort_interpro.r \
        -i $interpro_results \
        -u $eggnog_unmatched \
        -o ${params.pref}_interpro_matched.tsv \
        -f remaining_unmatched-${params.pref}.tsv
    """
    //
}
