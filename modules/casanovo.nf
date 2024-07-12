process CASANOVO {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"
    memory "20 GB"
    conda "/home/shannc/anaconda3/envs/casanovo"

    input:
    path(mzMLs)
    val(outdir)
    //
    output:
    path("${mzMLs.baseName}_casanovo.tsv"), emit: peps
    path("*.log"), emit: log
    //

    script:
    """
    casanovo \
        sequence \
        --config ${params.config_dir}/casanovo.yaml \
        --model $params.casanovomodel \
        --output temp.mztab \
        $mzMLs  > casanovo.log
    grep -v ^M temp.mztab > ${mzMLs.baseName}_casanovo.tsv
    """
    //
}

process EXTRACT_CASANOVO {
    publishDir "$outdir", mode: "copy"

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
    merge_tables.sh -r "$header" \
        -o all_casanovo.temp \
        -p tsv
    extract_denovo.py all_casanovo.temp casanovo
    """
    //
}
