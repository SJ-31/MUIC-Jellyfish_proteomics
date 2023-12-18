process INTERPROSCAN {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_fasta)
    val(outdir)
    //

    output:
    path("${unknown_fasta.baseName}-SCAN.tsv")
    //

    shell:
    def check = file("${outdir}/${unknown_fasta.baseName}-SCAN.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/!{unknown_fasta.baseName}-SCAN.tsv .
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
        cat header.txt temp.tsv > !{unknown_fasta.baseName}-SCAN.tsv
        '''
    }
    //
}

process SORT_INTERPRO {
    publishDir "$outdir", mode: "copy"

    input:
    path(interpro_results)
    path(eggnog_unmatched)
    val(outdir)
    //

    output:
    tuple path("interpro_anno-${name}.tsv"), path("interpro_meta-${name}.tsv")
    //

    script:
    name = interpro_results.baseName.replaceFirst(/-SCAN/, "")
    """
    Rscript $params.bin/sort_interpro.r \
        -i $interpro_results \
        -u $eggnog_unmatched \
        -m interpro_meta-${name}.tsv \
        -a interpro_anno-${name}.tsv \
    """
    //
}
