import pandas as pd
import sys

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import get_distances as gd

uniprot_downloads = (
    "/home/shannc/workflow/data/protein_databases/comparison_taxa"
)

sample_df = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"

sample_embd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/C_indra_esm_embeddings/embeddings.hdf5"
sample_embd_dist = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"

this = gd.hdf5ToDf(sample_embd)
df = pd.read_csv(sample_df, sep="\t")
filtered = df.query("category == 'venom_component'")
f = this[this.index.isin(filtered["ProteinId"])]

# saved = gd.getSaved(sample_embd, sample_embd_dist)
