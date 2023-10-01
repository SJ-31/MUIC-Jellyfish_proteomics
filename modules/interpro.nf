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
    '''
    header="query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways"
    interproscan.sh -i !{unknown_fasta} \
        -f tsv \
        -goterms \
        -pa \
        -b temp
    echo -e "$header" > header.txt
    cat header.txt temp.tsv > !{unknown_fasta.baseName}-SCAN
    '''
    //
}
