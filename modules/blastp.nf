process BLASTP {
    publishDir "$outdir", mode: "copy"

    input:
    path(query)
    val(db)
    val(outdir)
    //

    output:
    path("${query.baseName}-blast.csv")
    //

    shell:
    // Custom header format for blast tsv output
    //
    // "6 qseqid sseqid sstart send length bitscore evalue pident nident mismatch gaps"
    // qseqid, sseqid = query, subject seq id
    // sstart, send = start, end of alignment in subject
    // length = alignment length
    '''
    header="queryID,subjectID,sAlignStart,sAlignEnd,alignLen,bitscore,evalue,pident,nident,nmismatch,ngaps,"
    blastp -query !{query} \
        -db !{db} \
        -outfmt "10 delim=, qseqid sseqid sstart send length bitscore evalue pident nident mismatch gaps" | \
       sed "s;sp|\\(.*\\)|;\\1;" > temp.csv
    echo $header | cat - temp.csv > !{query.baseName}-blast.csv
    '''
    //
}
