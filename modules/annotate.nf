process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    conda "/home/shannc/anaconda3/envs/eggnog"
    // May need to install pandas and requests

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    tuple path("${name}_anno.tsv"), path("${name}_meta.tsv"), emit: annotations
    path("*fasta"), optional: true
    //

    shell:
    name = combined_tsv.baseName.replaceFirst( /db_hits-/, "")
    check = file("${outdir}/${name}_anno.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/!{name}_*.tsv .
        '''
    } else {
        template 'annotate.sh'
    }
    // need_anno = file("needs_annotating.fasta")
    // Some NCBI proteins fail to map to UniProt, so these are extracted,
    // annotated with InterProscan, then merged back into the original
    // annotations file
    // if (need_anno.exists()) {
    //     '''
    //     header='query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways'
    //     annotate.py -s \
    //         -u \
    //         -i annotated.tsv \
    //         -m !{name}_meta.tsv \
    //         -a !{name}_anno.tsv
    //     '''
        // need_anno = file("still_needs_annotating.fasta")
        // if (need_anno.exists()) {

        // }
}
