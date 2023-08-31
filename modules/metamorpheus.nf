process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "${params.pref}_results*"
    memory "10 GB"
    stageInMode "copy"
    // debug true
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(mzmls)
    val(outdir)
    val(database)
    //

    output:
    path("metamorpheus*")
    path("metamorpheus_AllPSMs.psmtsv"), emit: psms
    tuple val("metamorpheus"), path("metamorpheus_AllPSMs_FormattedForPercolator.tab"), emit: percolator
    //

    shell:
    '''
    metamorpheus -s !{mzmls} \
        -o . \
        -t !{params.config}/metamorpheus_params.toml \
        -d !{database}
    mv Task1SearchTask/[Aa]ll* .
    for i in [Aa]ll*
        do
        mv $i metamorpheus_${i}
        done
    mv metamorpheus_AllPSMs_FormattedForPercolator.tab edits.tab
    cat edits.tab | sed \
        -e 's/DECOY_/rev_/g' \
        -e 's/\[Common Variable:Oxidation on M\]/[15.9949]/g' \
        -e 's/\[Common Fixed:Carbamidomethyl on C\]//g' \
        > metamorpheus_AllPSMs_FormattedForPercolator.tab
    '''
}
