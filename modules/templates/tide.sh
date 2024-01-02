cp ${database} database.fasta

crux tide-index database.fasta db \
    --decoy-format peptide-reverse \
    --decoy-prefix rev_ \
    --nterm-protein-mods-spec 1X+42.010565 \
    --mods-spec 3M+15.994915

crux tide-search ${mgf} ./db \
    --auto-precursor-window warn \
    --spectrum-parser pwiz \
    --output-dir .

crux percolator ./tide-search.target.txt \
    --overwrite T \
    --decoy-prefix rev_ \
    --picked-protein database.fasta \
    --output-dir . \
    --protein-report-duplicates T \
    --search-input concatenated \
    --fileroot tide

rename_tide.sh
mv tide-search.target.txt tide_search.target.txt
mv tide-decoy.target.txt tide_search.decoy.txt
mv tide-search.log tide_search.log
