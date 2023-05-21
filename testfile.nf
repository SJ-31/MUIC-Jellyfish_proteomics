mqfiles = "resources/raw/CiCs*/combined/txt/allPeptides.txt"

test = Channel.fromPath(mqfiles).collectFile(name: "$projectDir/collect.txt", keepHeader: true, skip: 1)
// test = Channel.fromPath(mqfiles).collectFile(name: "$projectDir/collect.txt")