#!/usr/bin/env python
import sys

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/")

import map2go as mg

args = {
    "combined": "/home/shannc/2024-04-05_premap2go.tsv",
    "pfam2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam2go",
    "interpro2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/interpro2go",
    "ec2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/ec2go",
    "kegg2go": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/kegg_reaction2go",
    "pfam_db": "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam_entries.tsv",
}

combined = pd.read_csv(args["combined"], sep="\t")
before = combined.loc[:, ["GO"]].info()


# mapped = mapAllDb(
#     pfam_db_path=args["pfam_db"],
#     p2g_path=args["pfam2go"],
#     i2g_path=args["interpro2go"],
#     ec2g_path=args["ec2go"],
#     k2g_path=args["kegg2go"],
#     to_annotate=combined,
# )

# after = mapped.loc[:, ["GO"]].info()
