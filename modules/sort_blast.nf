process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(unknown_tsv), path(query_tsv) // Peptides known only from denovo peptides or from
    // transcriptome peptides
    path(database_hits) // Combined database proteins to merge accepted blast
    // results into
    path(blast_results)
    path(fasts_results)
    path(mapping) // Mapping file of ProteinId->Sequence, for retrieving metadata about
    // proteins newly inferred with blast
    val(outdir)
    //

    output:
    tuple path("${params.pref}_blast_unmatched.fasta"), \
        path("${params.pref}_blast_unmatched.tsv"), emit: unmatched
    path("${params.pref}_blast_matched.tsv"), emit: matched
    path("accepted_queries.tsv")
    //

    shell:
    check = file("${outdir}/${params.pref}_blast_matched.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/*_blast_* .
        cp !{outdir}/accepted_queries.tsv .
        '''
    } else {
        '''
        b_header="queryId,subjectId,sAlignStart,sAlignEnd,alignLen,bitscore,evalue,pident,nident,nmismatch,ngaps"
        f_header="queryId,subjectId,pident,alignLen,nmismatch,ngaps,qAlignStart,qAlignEnd,sAlignStart,aAlignEnd,evalue,bitscore"

        echo $b_header | cat - !{blast_results} > blast_results.csv
        echo $f_header | cat - !{fasts_results} > fasts_results.csv

        sort_blast.py -b blast_results.csv fasts_results.csv \
            --unknown_hits !{unknown_tsv} \
            --blast_query_map !{query_tsv} \
            -d !{database_hits} \
            -m !{mapping} \
            -i 85 \
            -e 0.00001 \
            -p 0.05 \
            -o !{params.pref}_blast_matched.tsv \
            -f !{params.pref}_blast_unmatched.fasta  \
            --unmatched_output !{params.pref}_blast_unmatched.tsv \
            --accepted_output accepted_queries.tsv
        '''
    // We want the tsv file for the unmatched peptides so that they their
    // metadata (which engines matched them, PEP, evalue etc.) will not be lost
    }
}

        // blastp header
        // header="queryId,subjectId,sAlignStart,sAlignEnd,alignLen,bitscore,evalue,pident,nident,nmismatch,ngaps"


        // fasts header
        // header="queryId,subjectId,pident,alignLen,nmismatch,ngaps,qAlignStart,qAlignEnd,sAlignStart,aAlignEnd,evalue,bitscore"
