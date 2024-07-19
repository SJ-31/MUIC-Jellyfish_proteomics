#!/usr/bin/env python
import sys

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/")

wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
ref = f"{wd}/data/reference"
tz = "/home/shannc/Downloads/thesis_testzone"
go_obo = f"{ref}/go.obo"
go_data = f"{ref}/go_data.tsv"
saved_go = f"{ref}/go_complete_saved.gml"

import map2go as mg
import go_subset as gs

args = {
    "combined": "/home/shannc/2024-04-05_premap2go.tsv",
    "pfam2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam2go",
    "interpro2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/interpro2go",
    "ec2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/ec2go",
    "kegg2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/kegg_reaction2go",
    "pfam_db": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam_entries.tsv",
}

# mapped = mapAllDb(
#     pfam_db_path=args["pfam_db"],
#     p2g_path=args["pfam2go"],
#     i2g_path=args["interpro2go"],
#     ec2g_path=args["ec2go"],
#     k2g_path=args["kegg2go"],
#     to_annotate=combined,
# )

# after = mapped.loc[:, ["GO"]].info()

GO = gs.CompleteGO(go_obo, go_data)


def test_successors1():
    success = gs.get_successors(
        GO.G,
        with_term=True,
        from_node="GO:0051179",
        ignore=["GO:0051640"],
        flatten=True,
    )
    assert "GO:0098953 receptor diffusion trapping" in success
    assert "GO:0050000 chromosome localization" not in success
    assert "GO:0051640 organelle localization" not in success
    assert "GO:0051653 spindle localization" not in success


def test_successors2():
    success = gs.get_successors(
        GO.G,
        with_term=True,
        from_node="GO:0032774",
        ignore=["GO:0001172", "GO:0039694"],
        flatten=True,
    )
    assert "GO:0001172 RNA-templated transcription" not in success
    assert "GO:0140745 siRNA transcription" not in success
    assert "GO:0006362 transcription elongation by RNA polymerase I" in success
    assert (
        "GO:0001175 transcriptional start site selection at RNA polymerase III promoter"
        in success
    )
