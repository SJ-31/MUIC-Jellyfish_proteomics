process TIDE {
    publishDir "$outdir", mode: "copy", pattern: "tide-search*"
    publishDir "$percolatordir", mode: "copy", pattern: "*percolator*"

    input:
    path(mzXMLs)
    val(outdir)
    val(percolatordir)
    //

    output:
    path("tide-search.target.txt")
    path("tide-search.decoy.txt")
    path("*percolator*"), emit: percolator
    // ADd a special output channel for percolator

    script:
    """
    crux tide-index $params.database db \
         --decoy-format peptide-reverse \
         --decoy-prefix rev_ \
        --nterm-protein-mods-spec 1X+42.010565 \
        --mods-spec 3M+15.994915

    crux tide-search $mxXMLs ./db \
        --auto-precursor-window warn \
        --spectrum-parser mstoolkit \
        --output-dir .

    crux percolator ./tide-search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --picked-protein $params.database \
    --output-dir . \
    --protein-report-duplicates T \
    --search-input concatenated \
    --fileroot tide
    """
    //
}
