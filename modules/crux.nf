process TIDE {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzXMLs)
    //

    output:

    //

    script:
    """
    crux tide-index $params.database db \
         --decoy-format peptide-reverse \
         --decoy-prefix rev_

    crux tide-search $mxXMLs ./db \
        --auto-precursor-window warn \
        --spectrum-parser mstoolkit \

    crux percolator ./tide-search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --

    """
    //
}
