process MSFRAGGER {
    memory { 2.GB * task.attempt }
    publishDir "$outdir", mode: "copy"

    input:
    path(mzmls)
    val(pars)
    val(outdir)

    output:
    path("${params.pref}_msfragger.tsv")
    path("msfragger*.tsv")
    path("${params.pref}-msfragger_percolator")

    script:
    tsv_header = "scannum	precursor_neutral_mass	retention_time	charge	ion_mobility	compensation_voltage	hit_rank	peptide	peptide_prev_aa	peptide_next_aa	protein	num_matched_ions	tot_num_ions	calc_neutral_pep_mass	massdiff	num_tol_term	num_missed_cleavages	modification_info	hyperscore	nextscore	expectscore	best_locs	score_without_delta_mass	best_score_with_delta_mass	second_best_score_with_delta_mass	delta_score"
    pin_header = "SpecId	Label	ScanNr	ExpMass	retentiontime	rank	abs_ppm	isotope_errors	log10_evalue	hyperscore	delta_hyperscore	matched_ion_num	complementary_ions	ion_series	weighted_average_abs_fragment_ppm	peptide_length	ntt	nmc	charge_1	charge_2	charge_3	charge_4	charge_5	charge_6	charge_7_or_more	15.994915M	Peptide	Proteins"
    """
    java -Xmx32g -jar ~/tools/MSFragger-3.7/MSFragger-3.7.jar  \
        ${pars} \
        ${mzmls}

    merge_tables.sh -r "$tsv_header" \
        -o ${params.pref}_msfragger.txt \
        -p tsv
    mv ${params.pref}_msfragger.txt ${params.pref}_msfragger.tsv

    percolator_wrapper.sh \
       -p $params.pref \
       -i pin \
       -f $params.database \
       -e msfragger
    """
    // The output pin files need to be merged into a single one before doing the percolator analysis
    // Calibrate mass
}

