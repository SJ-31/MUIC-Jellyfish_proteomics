process MAXQUANT {
    stageInMode "copy"
    publishDir "$outdir", mode: "copy"

    input:
    path(raw_file)
    val(mq_config)
    val(outdir)
    //
    output:
    path("${raw_file.baseName}_combined")

    shell:
    '''
    cp !{mq_config} template.xml
    maxquant_wrapper.py !{raw_file} template.xml
    rm template.xml
    dotnet /home/shannc/tools/MaxQuant_2.4.2.0/bin/MaxQuantCmd.exe mqconfig.xml
    mv combined/txt ./!{raw_file.baseName}_combined
    '''
    //
}
