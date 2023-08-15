process TIDE {
    publishDir "$outdir", mode: "copy", pattern: "tide-search*"
    publishDir "$percolatordir", mode: "copy", pattern: "*percolator*"

    input:
    path(mzXMLs)
    val(outdir)
    val(percolatordir)
    val(database)
    //

    output:
    path("tide-search.target.txt")
    path("tide-search.decoy.txt")
    path("*percolator*"), emit: percolator
    // ADd a special output channel for percolator

    script:
    """
    cp $database database.fasta
    crux tide-index database.fasta db \
         --decoy-format peptide-reverse \
         --decoy-prefix rev_ \
        --nterm-protein-mods-spec 1X+42.010565 \
        --mods-spec 3M+15.994915

    crux tide-search $mzXMLs ./db \
        --auto-precursor-window warn \
        --spectrum-parser mstoolkit \
        --output-dir .

    crux percolator ./tide-search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --picked-protein database.fasta \
    --output-dir . \
    --protein-report-duplicates T \
    --search-input concatenated \
    --fileroot tide
    """
    //
}

process TIDE_COMBINED_PEP {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_files)
    val(outdir)
    //

    output:
    path("tide_psm2combined_PEP.tsv"), emit: psm2combinedPEP
    path("tide_prot2combined_PEP.tsv"), emit: prot2combinedPEP
    //

    script:
    """
    Rscript $params.bin/to_combined_PEP.r \
        -m tide.percolator.target.psms.txt \
        -d tide.percolator.decoy.psms.txt \
        -o tide_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m tide.percolator.target.proteins.txt \
        -d tide.percolator.decoy.proteins.txt \
        --protein_matches \
        -o tide_prot2combined_PEP.tsv
    """
    //
}
