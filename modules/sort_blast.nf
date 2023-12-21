process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_tsv) // Peptides known only from denovo peptides or from
    // transcriptome peptides
    path(unmatched_peptides_tsv) // Peptides not matched during database search
    path(database_hits) // Combined database proteins to merge accepted blast
    // results into
    path(blast_results)
    path(mapping)
    path(blast_query)
    each arg_string
    val(outdir)
    //

    output:
    tuple path("${prefix}_unmatched.fasta"), path("${prefix}_unmatched.tsv"), emit: unmatched
    path("db_hits-${prefix}.tsv"), emit: matched
    //

    script:
    // Syntax is <prefix> <do_one_hit?> <do_best_only?> <identity_threshold> <pep_threshold> <evalue_threshold>
    arguments = arg_string.split(/ /)
    prefix = arguments[0]
    one_hit = arguments[1] // 0 Removes proteins identified by only one peptide
    best = arguments[2] // 1 Keeps best hit only, no degenerate peptides allowed
    identity_thresh = arguments[3]
    p_thresh = arguments[4]
    evalue_thresh = arguments[5]
    """
    grep ">" ${blast_query} | sed 's/>//' > blast_query.txt
    sort_blast.py -b $blast_results \
        --unknown_hits $unknown_tsv \
        --unmatched_peptides $unmatched_peptides_tsv \
        -d $database_hits \
        -m $mapping \
        -q blast_query.txt \
        -f "${prefix}_unmatched.fasta"  \
        --unmatched_tsv ${prefix}_unmatched.tsv \
        -i $identity_thresh \
        -e $evalue_thresh \
        -p $p_thresh \
        --one_hit $one_hit \
        --keep_best $best \
        -o db_hits-${prefix}.tsv
    """
    // We want the tsv file for the unmatched peptides so that they their
    // metadata (which engines matched them, PEP, evalue etc.) will not be lost
}
