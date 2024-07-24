process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_tsv) // Peptides known only from denovo peptides or from
    // transcriptome peptides
    path(database_hits) // Combined database proteins to merge accepted blast
    // results into
    path(blast_results)
    path(mapping) // Mapping file of ProteinId->Sequence, for retrieving metadata about
    // proteins newly inferred with blast
    path(blast_query)
    val(outdir)
    //

    output:
    tuple path("${params.pref}_blast_unmatched.fasta"), \
        path("${params.pref}_blast_unmatched.tsv"), emit: unmatched
    path("${params.pref}_blast_matched.tsv"), emit: matched
    //

    shell:
    check = file("${outdir}/${params.pref}_blast_matched.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/*_blast_* .
        '''
    } else {
        '''
        header="queryID,subjectID,pident,alignLen,nmismatch,ngaps,qAlignStart,qAlignEnd,sAlignStart,aAlignEnd,evalue,bitscore,"

        echo $header | cat - !{blast_results} > blast_results.csv
        grep ">" !{blast_query} | sed 's/>//' > blast_query.txt

        sort_blast.py -b blast_results.csv \
            --unknown_hits !{unknown_tsv} \
            -d !{database_hits} \
            -m !{mapping} \
            -q blast_query.txt \
            -f !{params.pref}_blast_unmatched.fasta  \
            --unmatched_tsv !{params.pref}_blast_unmatched.tsv \
            -i 80 \
            -e 0.00001 \
            -p 0.05 \
            -o !{params.pref}_blast_matched.tsv
        '''
    // We want the tsv file for the unmatched peptides so that they their
    // metadata (which engines matched them, PEP, evalue etc.) will not be lost
    }
}

        // OLD unused blast header
        // header="queryID,subjectID,sAlignStart,sAlignEnd,alignLen,bitscore,evalue,pident,nident,nmismatch,ngaps,"
