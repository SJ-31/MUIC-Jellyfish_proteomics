process BLASTP {
    publishDir "$outdir", mode: "copy"

    input:
    path(query)
    val(db)
    tuple val(min), val(max)
    val(outdir)
    //

    output:
    path(blast_results)
    //

    shell:
    // Custom header format for blast tsv output
    //
    // "6 qseqid sseqid sstart send length bitscore evalue pident nident mismatch gaps"
    // qseqid, sseqid = query, subject seq id
    // sstart, send = start, end of alignment in subject
    // length = alignment length
    blast_results = "${query.baseName}${min}-blast.csv"
    check = file("${outdir}/${blast_results}")
    if (check.exists()) {
        '''
        cp !{outdir}/!{blast_results} .
        '''
    } else {
        if (max == 35) {
            flags = "-task blastp-short"
        } else if (max == 50) {
            flags = "-task blastp -matrix PAM70"
        } else if (max == 85) {
            flags = "-task blastp -matrix BLOSUM80"
        } else {
            flags = "-task blastp -matrix BLOSUM62"
        }
        '''
        seqkit seq -m !{min} -M !{max} !{query} > filtered.fasta
        blastp -query filtered.fasta \
            !{flags} \
            -db !{db} \
            -outfmt "10 delim=, qseqid sseqid sstart send length bitscore evalue pident nident mismatch gaps" | \
            sed "s;sp|\\(.*\\)|;\\1;" > !{blast_results}
        '''
    }
    //
}

