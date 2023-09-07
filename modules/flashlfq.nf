process FLASHLFQ {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "*{txt,toml}"
    errorStrategy 'ignore'

    input:
    path(scan_prot_mappings)
    path(mzmls)
    val(outdir)
    //

    output:
    path("Quantified*")
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
