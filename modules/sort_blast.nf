process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_hits) // Peptides known only from denovo peptides or from
    // transcriptome peptides
    path(database_hits)
    path(blast_results)
    path(mapping)
    path(unmatched_peptides) // Peptides not matched during database search
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
    sort_blast.py -b $blast_results \
        --unknown_hits $unknown_hits \
        -d $database_hits \
        -m $mapping \
        -f unmatched.fasta  \
        --unmatched_tsv ${prefix}_unmatched.tsv \
        -i $identity_thresh \
        -e $evalue_thresh \
        -p $p_thresh \
        --one_hit $one_hit \
        --keep_best $best \
        -o db_hits-${prefix}.tsv
    cat $unmatched_peptides unmatched.fasta > temp.fasta
    cd-hit -i temp.fasta -o ${prefix}_unmatched.fasta -c 1.0
    """
    // We want the tsv file for the unmatched peptides so that they their
    // metadata (which engines matched them, PEP, evalue etc.) will not be lost
}
