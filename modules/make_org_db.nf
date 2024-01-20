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
    splits = params.species_spec.split("_")
    name = annotations.baseName.replaceFirst(/_meta/, "")
    """
    Rscript $params.bin/make_orgdb.r \
        -c $annotations \
        -m "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -a "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -t $params.tax_id \
        -s ${splits[0]} \
        -g ${splits[1]}

    mv org* ${name}_org.db
    """
    //
}
