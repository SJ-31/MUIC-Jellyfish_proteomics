process FLASHLFQ {
    publishDir "$outdir", mode: "copy", pattern: "*tsv"
    publishDir "$logdir", mode: "copy", pattern: "*{txt,toml}"

    input:
    path(scan_prot_mappings)
    path(mzmls)
    val(logdir)
    val(outdir)
    //

    output:
    path("Quantified*")
    path("QuantifiedProteins.tsv"), emit: prot
    path("*.{txt,toml}")
    path("flashlfq.tsv")
    //

    script:
    """
    mkdir scan_path; mv $scan_prot_mappings scan_path
    mkdir mzML; mv *.mzML mzML

    Rscript $params.bin/R/flashlfq_format2.r \
        -p scan_path \
        -o flashlfq.tsv


    $params.dotnet6 $params.flashlfq \
        --idt flashlfq.tsv \
        --out . \
        --rep mzML \
        --ppm 5
    """
    // TODO: Put the "ExperimentalDesign.tsv" file into 'scan_path' after creating it
    // with the Rscript
}
