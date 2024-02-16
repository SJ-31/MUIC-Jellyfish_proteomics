process MAKE_ORGDB {
    publishDir "$outdir", mode: "copy"

    input:
    path(annotations)
    val(comparison_data)
    val(outdir)
    //

    output:
    path("${name}_org.db")
    //

    script:
    splits = params.species_spec.split("_")
    name = annotations.baseName.replaceFirst(/_meta/, "")
    """
    Rscript $params.bin/R/make_orgdb.r \
        -c $annotations \
        -m "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -a "Shann Chongwattananukul <shann.cho@student.mahidol.edu>" \
        -p $comparison_data \
        -t $params.tax_id \
        -s ${splits[0]} \
        -g ${splits[1]}

    mv org* ${name}_org.db
    """
    //
}
