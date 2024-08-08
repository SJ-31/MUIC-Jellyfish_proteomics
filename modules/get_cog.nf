process GET_COG {
    publishDir "$outdir", mode: "copy"
    input:
    path(combined_tsv)
    val(outdir)
    //
    output:
    path(output), emit: tsv
    path("cog_assignment_evidence.tsv")
    //
    script:
    output = combined_tsv.baseName.replace("wcoverage", "wcog.tsv")
    """
    grouping.py -i $combined_tsv \
        -o $output \
        --pfam_map_path $params.pfam_clans \
        --cog_map_path $params.cog_maps \
        --all_map_path $params.all_maps \
        --toxin_map_path $params.toxin_maps \
        --evidence_file_output cog_assignment_evidence.tsv \
        --mode cog
    """
    //
}
