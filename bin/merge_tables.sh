#!/usr/bin/env sh

while getopts "r:h:o:p:" opt; do
    case $opt in
        r)
            header=$OPTARG ;;
        o)
            output=$OPTARG ;;
        p)
            pattern=$OPTARG ;;
        h)
            echo "-r <header> -o <output_filename> -p <pattern_of_files_to_merge>"
    esac
done


echo "$header" > "$output"
find . -maxdepth 1 -name "*${pattern}" -exec sed -se 1d {} + >> "$output"
