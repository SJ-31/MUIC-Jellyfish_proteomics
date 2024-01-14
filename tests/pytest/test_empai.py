#!/usr/bin/env ipython

import sys

sys.path.append("./bin")
import emPAI as ep

pth = "./results/jellyfish/1-First_pass/"
test = f"{pth}/Unmatched/Database-annotated/jellyfish_downloads_anno-3.tsv"
df = pd.read_csv(test, sep="\t")
