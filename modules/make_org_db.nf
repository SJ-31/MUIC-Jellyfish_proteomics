process MAKE_ORGDB {
    publishDir "$outdir", mode: "copy"

    input:
    path(annotations)
    val(outdir)
    //

    output:
    path("${name}_org.db")
    //

    script:
    name = annotations.baseName.replaceFirst(/_meta/, "")
    """
    Rscript $params.bin/make_orgdb.r \
        -c $annotations \
        -m "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -a "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -t $params.tax_id \
        -s $params.species \
        -g $params.genus

    mv org* ${name}_org.db
    """
    //
}
