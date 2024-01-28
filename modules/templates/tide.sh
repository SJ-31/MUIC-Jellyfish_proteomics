cp !{database} database.fasta

tsv_header="file	scan	charge	spectrum precursor m/z	spectrum neutral mass	peptide mass	delta_cn	delta_lcn	xcorr score	xcorr rank	distinct matches/spectrum	sequence	modifications	cleavage type	protein id	flanking aa	target/decoy"

crux tide-index database.fasta db \
    --decoy-format peptide-reverse \
    --decoy-prefix rev_ \
    --clip-nterm-methionine T \
    --keep-terminal-aminos none \
    --min-length 7 \
    --max-length 50 \
    --nterm-protein-mods-spec 1X+42.010565 \
    --mods-spec 3M+15.994915 \
    --enzyme trypsin 
    # --missed-cleavages 2

for mgf_file in *mgf
do
    name=$(echo "$mgf_file" | sed 's/.mgf//')
    mkdir "${name}"; mv "$mgf_file" "${name}"; cd "${name}"
    crux tide-search "${mgf_file}" ../db \
        --spectrum-parser pwiz \
        --output-dir .
    mv tide-search.target.txt ../"${name}".target.txt
    mv tide-search.decoy.txt ../"${name}".decoy.txt
    mv tide-search.log.txt ../tide_search_"${name}".log
    cd ..
done


merge_tables.sh -r "$tsv_header" \
    -o tide_search.target.temp \
    -p target.txt
mv tide_search.target.temp tide_search.target.txt # Will cause percolator error otherwise

merge_tables.sh -r "$tsv_header" \
    -o tide_search.decoy.temp \
    -p decoy.txt
mv tide_search.decoy.temp tide_search.decoy.txt

crux percolator ./tide_search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --picked-protein database.fasta \
    --output-dir . \
    --protein-report-duplicates T \
    --search-input concatenated \
    --fileroot tide

for t in tide*percolator*
    do
    rename="$(echo "$t" | sed -e 's/\./_/g' -e 's/_txt/\.tsv/' -e 's/_target//')"
    mv "$t" "$rename"
done

# rename_tide.sh
