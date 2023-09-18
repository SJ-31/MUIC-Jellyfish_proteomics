process SORT_BLAST {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_hits)
    path(database_hits)
    path(blast_results)
    path(mapping)
    val(one_hit) // 0 Removes proteins identified by only one peptide
    val(best) // 1 Keeps best hit only, no degenerate peptides allowed
    val(identity_thresh)
    val(pep_thresh)
    val(evalue_thresh)
    val(prefix)
    val(outdir)
    //

    output:
    path("${prefix}_unmatched.fasta")
    path("blast_matched-${prefix}.tsv")
    //

    script:
    """
    sort_blast.py -b $blast_results \
        -u $unknown_hits \
        -d $database_hits \
        -m $mapping \
        -f ${prefix}_unmatched.fasta \
        -i $identity_thresh \
        -e $evalue_thresh \
        -p $pep_thresh \
        --one_hit $one_hit \
        --keep_best $best \
        -o blast_matched-${prefix}.tsv
    """
    //
}
