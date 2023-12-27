process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    publishDir "$outdir", mode: "copy", pattern: "${name}_eggnog-unmatched"
    publishDir "$outdir", mode: "copy", pattern: "${name}_interpro-unmatched"
    // May need to install pandas and requests

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    tuple path("${name}_anno.tsv"), path("${name}_meta.tsv"), emit: annotations
    path("${name}_eggnog-unmatched")
    path("${name}_interpro-unmatched")
    path("*fasta"), optional: true
    //

    shell:
    name = combined_tsv.baseName.replaceFirst( /db_hits-/, "")
    template 'annotate.sh'
}
