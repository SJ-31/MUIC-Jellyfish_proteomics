process FLASHLFQ {
    publishDir "$outdir", mode: "copy"
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
    Rscript $params.bin/flashlfq_format2.r \
        -p scan_path \
        -o flashlfq.tsv

    $params.dotnet6 $params.flashlfq \
        --idt flashlfq.tsv \
        --out . \
        --rep mzML \
        --ppm 5
    """
    //
}
