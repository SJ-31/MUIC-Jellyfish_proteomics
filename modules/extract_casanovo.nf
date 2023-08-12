process EXTRACT_CASANOVO {
    publishDir "$output", mode: "copy"

    input:
    path(casanovo_output)
    val(outdir)
    //

    output:
    path("casanovo*fasta")
    //

    script:
    header = "PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end	opt_ms_run[1]_aa_scores"
    """
    merge_tables.sh -r $header \
        -o all_casanovo.temp \
        -p tsv
    extract_casanovo.py all_casanovo.temp
    """
    //
}
