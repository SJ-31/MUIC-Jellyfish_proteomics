process MAXQUANT {

    stageInMode "copy"
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(raw_files) // Maxquant seems to only work with .raw files
    val(mq_config)
    val(outdir)
    val(logdir)
    val(database)
    //
    output:
    path("all_combined")
    path("msms.txt"), emit: msms
    path("msmsScans.txt"), emit: msmsScans
    path("*.log")

    shell:
    '''
    cp !{params.config}/!{mq_config} template.xml
    maxquant_wrapper.py -c template.xml \
        -r . -o mqconfig.xml \
        -d !{database} \
        --ptms False --fraction 32767 --index 0 \
        --experiment "" --reference_channel ""
    rm template.xml
    dotnet /home/shannc/tools/MaxQuant_2.4.2.0/bin/MaxQuantCmd.exe mqconfig.xml \
        > maxquant.log
    mv combined/txt ./all_combined
    cp all_combined/msms.txt msms.txt
    cp all_combined/msmsScans.txt msmsScans.txt
    rm -R combined/search
    '''
    //
    //
    // cp !{params.config}/maxquant_ms2rescore.xml template.xml Use this for ms2rescore
}


// process MAXQUANT {

//     stageInMode "copy"
//     publishDir "$outdir", mode: "copy"
//     publishDir "$params.logs", mode: "copy", pattern: "*.log"

//     input:
//     path(raw_file) // Maxquant seems to only work with .raw files
//     val(outdir)
//     val(database)
//     //
//     output:
//     path("${raw_file.baseName}_combined")
//     path("${raw_file.baseName}_msms.txt"), emit: msms
//     path("${raw_file.baseName}_msmsScans.txt"), emit: msmsScans
//     path("*.log")

//     shell:
//     '''
//     cp !{params.config}/maxquant_ms2rescore.xml template.xml
//     maxquant_wrapper.py !{raw_file} template.xml !{database}
//     rm template.xml
//     dotnet /home/shannc/tools/MaxQuant_2.4.2.0/bin/MaxQuantCmd.exe mqconfig.xml \
//         > maxquant_!{raw_file.baseName}.log
//     mv combined/txt ./!{raw_file.baseName}_combined
//     cp !{raw_file.baseName}_combined/msms.txt !{raw_file.baseName}_msms.txt
//     cp !{raw_file.baseName}_combined/msmsScans.txt !{raw_file.baseName}_msmsScans.txt
//     '''
//     //
//     //
//     // cp !{params.config}/maxquant_ms2rescore.xml template.xml Use this for ms2rescore
// }
