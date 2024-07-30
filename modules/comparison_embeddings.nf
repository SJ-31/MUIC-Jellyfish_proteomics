process COMPARISON_EMBEDDINGS {
    publishDir "$outdir", mode: "copy"

    input:
    path(sample_embeddings)
    path(sample_tsv)
    path(comparison_embeddings)
    path(comparison_tsv)
    val(embedding_type)
    val(outdir)
    //

    output:
    path("sample_with_comp_dist.hdf5"), emit: dist
    path("sample_with_comp_embd.hdf5"), emit: embd
    //

    script:
    def check = file("$outdir/sample_with_comp_dist.hdf5")
    if (check.exists()) {
        """
        cp $outdir/*hdf5 .
        """
    } else if (embedding_type == "protein") {
        """
        get_distances.py \
            -i $sample_embeddings \
            -o sample_with_comp_dist.hdf5 \
            --sample_tsv $sample_tsv \
            --filter_criteria "category == 'venom_component'" \
            --comparison_embd $comparison_embeddings \
            --comparison_tsv $comparison_tsv \
            --write_embd_file sample_with_comp_embd.hdf5
        conda deactivate
        """
    } else if (embedding_type == "GO") {
        """
        Rscript $params.bin/R/analysis/prepare_embeddings.r \
            --r_source "$params.bin/R" \
            --python_source "$params.bin" \
            --sample_tsv $sample_tsv \
            --embedding_path $sample_embeddings \
            --sample_name $params.pref \
            --comparison_tsv $comparison_tsv \
            --comparison_embd $comparison_embeddings \
            --embd_output sample_with_comp_embd.hdf5 \
            --dist_output sample_with_comp_dist.hdf5
        """
    }
    //
}
