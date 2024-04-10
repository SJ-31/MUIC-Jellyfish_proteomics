process COMPARISON_EMBEDDINGS {
    conda "/home/shannc/anaconda3/envs/reticulate"
    publishDir "$outdir", mode: "copy"

    input:
    path(sample_embeddings)
    path(sample_tsv)
    path(comparison_embeddings)
    path(comparison_tsv)
    val(outdir)
    //

    output:
    path("sample_with_comp_dist.hdf5")
    path("sample_with_comp_embd.hdf5")
    //

    script:
    def check = file("$outdir/sample_with_comp_dist.hdf5")
    if (check.exists()) {
        """
        cp $outdir/*hdf5 .
        """
    } else {
        """
        get_distances.py -i $sample_embeddings \
            -o sample_with_comp_dist.hdf5 \
            --sample_tsv $sample_tsv \
            --filter_criteria "category == 'venom_component'" \
            --comparison_embd $comparison_embeddings \
            --comparison_tsv $comparison_tsv \
            --write_embd_file sample_with_comp_embd.hdf5
        """
    }
    //
}
