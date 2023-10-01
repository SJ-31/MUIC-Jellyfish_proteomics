process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_hits)
    path(database_hits)
    path(blast_results)
    path(mapping)
    path(unmatched_peptides)
    each arg_string
    val(outdir)
    //

    output:
    path("${prefix}_unmatched.fasta"), emit: unmatched
    path("blast_matched-${prefix}.tsv"), emit: matched
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
    sort_blast.py -b $blast_results \
        -u $unknown_hits \
        -d $database_hits \
        -m $mapping \
        -f unmatched.fasta  \
        -t ${prefix}_unmatched.tsv \
        -i $identity_thresh \
        -e $evalue_thresh \
        -p $p_thresh \
        --one_hit $one_hit \
        --keep_best $best \
        -o blast_matched-${prefix}.tsv
    cat $unmatched_peptides unmatched.fasta > ${prefix}_unmatched.fasta
    """
    //
}
