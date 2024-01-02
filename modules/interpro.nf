process INTERPROSCAN {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_fasta)
    val(outdir)
    //

    output:
    path("${name}_interpro.tsv")
    //

    shell:
    name = unknown_fasta.baseName.replaceFirst(/_eggnog/, "")
    check = file("${outdir}/${name}_interpro.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/!{name}_interpro.tsv .
        '''
    } else {
        '''
        header="query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways"
        interproscan.sh -i !{unknown_fasta} \
            -f tsv \
            -goterms \
            -pa \
            -b temp
        echo -e "$header" > header.txt
        cat header.txt temp.tsv > !{name}_interpro.tsv
        '''
    }
    //
}

process SORT_INTERPRO {
    publishDir "$outdir", mode: "copy", pattern: "*{meta,anno}*"
    publishDir "$outdir_unmatched", mode: "copy", pattern: "remaining_unmatched*"

    input:
    path(interpro_results)
    path(eggnog_unmatched)
    val(outdir)
    val(outdir_unmatched)
    //

    output:
    tuple val(name), path("interpro_anno-${name}.tsv"), path("interpro_meta-${name}.tsv"), emit: matched
    path("remaining_unmatched-${name}.tsv")
    //

    script:
    name = interpro_results.baseName.replaceFirst(/_interpro/, "")
    """
    Rscript $params.bin/sort_interpro.r \
        -i $interpro_results \
        -u $eggnog_unmatched \
        -m interpro_meta-${name}.tsv \
        -a interpro_anno-${name}.tsv \
        -f remaining_unmatched-${name}.tsv \
    """
    //
}
