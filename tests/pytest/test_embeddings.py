import sys

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import get_distances as gd

sample_embd = (
    "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/C_indra_esm_embeddings/embeddings.hdf5",
)
sample_embd_dist = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"

this = gd.getEmbeddings(sample_embd)

saved = gd.getSaved(sample_embd, sample_embd_dist)
