process TIDE {
    publishDir "$outdir", mode: "copy", pattern: "tide-search*"
    publishDir "$percolatordir", mode: "copy", pattern: "*percolator*"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    path(mgf)
    val(outdir)
    val(percolatordir)
    val(database)
    //

    output:
    path("tide-search.target.txt"), emit: target
    path("tide-search.decoy.txt"), emit: decoy
    path("*percolator*"), emit: percolator
    path("tide_percolator_proteins.tsv"), emit: perc_protein
    path("tide_percolator_psms.tsv"), emit: perc_psms
    path("*.log")
    // ADd a special output channel for percolator

    script:
    """
    cp ${database} database.fasta
    crux tide-index database.fasta db \
         --decoy-format peptide-reverse \
         --decoy-prefix rev_ \
        --nterm-protein-mods-spec 1X+42.010565 \
        --mods-spec 3M+15.994915

    crux tide-search ${mgf} ./db \
        --auto-precursor-window warn \
        --spectrum-parser mstoolkit \
        --output-dir . > tide.log

    crux percolator ./tide-search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --picked-protein database.fasta \
    --output-dir . \
    --protein-report-duplicates T \
    --search-input concatenated \
    --fileroot tide > percolator_tide.log

    rename_tide.sh
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
        -m tide_percolator_psms.tsv \
        -d tide_percolator_decoy_psms.tsv \
        -o tide_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m tide_percolator_proteins.tsv \
        -d tide_percolator_decoy_proteins.tsv \
        --protein_matches \
        -o tide_prot2combined_PEP.tsv
    """
    //
}
