process MAXQUANT {
    stageInMode "copy"
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    path(raw_file)
    val(outdir)
    val(database)
    //
    output:
    path("${raw_file.baseName}_combined")
    tuple val("maxquant"), path("${raw_file.baseName}_msms.txt"), emit: ms2rescore

    shell:
    '''
    cp !{params.config}/maxquant_ms2rescore.xml template.xml
    maxquant_wrapper.py !{raw_file} template.xml !{database}
    rm template.xml
    dotnet /home/shannc/tools/MaxQuant_2.4.2.0/bin/MaxQuantCmd.exe mqconfig.xml \
        > maxquant_!{raw_file.baseName}.log
    mv combined/txt ./!{raw_file.baseName}_combined
    cp !{raw_file.baseName}_combined/msms.txt !{raw_file.baseName}_msms.txt
    '''
    //
}
