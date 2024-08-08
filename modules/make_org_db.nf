process MAKE_ORGDB {
    publishDir "$outdir", mode: "copy"

    input:
    path(annotations)
    val(outdir)
    //

    output:
    path("${params.pref}_org.db")
    //

    script:
    splits = params.species_spec.split("_")
    """
    Rscript $params.bin/R/make_orgdb.r \
        -c $annotations \
        -m "$params.author" \
        -a "$params.author" \
        -t $params.tax_id \
        -s ${splits[0]} \
        -g ${splits[1]}

    mv org* ${params.pref}_org.db
    """
    //
}

