process MS2RESCORE {
    publishDir "$outdir", mode: "copy"
    memory "10 GB"
    stageInMode "copy"
    conda "/home/shannc/anaconda3/envs/ms2rescore"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    tuple val(engine), path(engine_out)
    val(outdir)
    path(mgfs)
    //

    output:
    tuple val(engine), path("${engine}_all_pins.temp"), emit: pin
    path("*.log")
    //

    shell:
    '''
    mkdir mgf
    mv *.mgf mgf
    pin_header="SpecId	Label	ScanNr	spec_pearson_norm	ionb_pearson_norm	iony_pearson_norm	spec_mse_norm	ionb_mse_norm	iony_mse_norm	min_abs_diff_norm	max_abs_diff_norm	abs_diff_Q1_norm	abs_diff_Q2_norm	abs_diff_Q3_norm	mean_abs_diff_norm	std_abs_diff_norm	ionb_min_abs_diff_norm	ionb_max_abs_diff_norm	ionb_abs_diff_Q1_norm	ionb_abs_diff_Q2_norm	ionb_abs_diff_Q3_norm	ionb_mean_abs_diff_norm	ionb_std_abs_diff_norm	iony_min_abs_diff_norm	iony_max_abs_diff_norm	iony_abs_diff_Q1_norm	iony_abs_diff_Q2_norm	iony_abs_diff_Q3_norm	iony_mean_abs_diff_norm	iony_std_abs_diff_norm	dotprod_norm	dotprod_ionb_norm	dotprod_iony_norm	cos_norm	cos_ionb_norm	cos_iony_norm	spec_pearson	ionb_pearson	iony_pearson	spec_spearman	ionb_spearman	iony_spearman	spec_mse	ionb_mse	iony_mse	min_abs_diff_iontype	max_abs_diff_iontype	min_abs_diff	max_abs_diff	abs_diff_Q1	abs_diff_Q2	abs_diff_Q3	mean_abs_diff	std_abs_diff	ionb_min_abs_diff	ionb_max_abs_diff	ionb_abs_diff_Q1	ionb_abs_diff_Q2	ionb_abs_diff_Q3	ionb_mean_abs_diff	ionb_std_abs_diff	iony_min_abs_diff	iony_max_abs_diff	iony_abs_diff_Q1	iony_abs_diff_Q2	iony_abs_diff_Q3	iony_mean_abs_diff	iony_std_abs_diff	dotprod	dotprod_ionb	dotprod_iony	cos	cos_ionb	cos_iony	observed_retention_time	predicted_retention_time	rt_diff	rt_diff_best	observed_retention_time_best	predicted_retention_time_best	RawScore	RawDeltaScore	RawModLocProb	ChargeN	Mass	PepLen	dM	enzInt	absdM	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	Charge7	MeanErrorTop7	sqMeanErrorTop7	StdevErrorTop7	lnExplainedIonCurrent	lnNTermIonCurrentRatio	lnCTermIonCurrentRatio	lnMS2IonCurrent	Peptide	Proteins"
    for output in !{engine_out}
       do
        name=$(echo $output | sed 's/_.*//')
        ms2rescore $output \
            -c !{params.config}/ms2rescore_config.json \
            -o $name \
            -m mgf >> ms2rescore.log
        done
    merge_tables.sh -r "${pin_header}" \
        -o !{engine}_all_pins.temp \
        -p pin

    case !{engine} in
        maxquant)
            add_flanks.py !{engine}_all_pins.temp ;;
        *)
            echo default ;;
    esac
    '''
    //
}
