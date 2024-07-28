process FASTS {
    publishDir "$outdir", mode: "copy"

    input:
    path(query)
    val(db)
    val(outdir)
    //

    output:
    path(fasts_results)
    //

    shell:
    fasts_results = "${query}".replace("unknown", "fasts").replace("fasta", "csv")
    check = file("${outdir}/${fasts_results}")
    if (check.exists()) {
        '''
        cp !{outdir}/!{fasts_results} .
        '''
    } else {
        '''
        fasts36 -m 8 !{query} !{db} \
            -E 0.001 | sed 's/\t/,/g' > !{fasts_results}
        '''
    }
    //
}
