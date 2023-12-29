#!/usr/bin/env ipython
#
import pandas as pd
import sys

pth = "./results/jellyfish/1-First_pass"
sys.path.append("./bin")
import merge_quantifications as mq

args = {
    "direct_lfq": f"{pth}/Quantify/directlfq_prot.tsv",
    "flash_lfq": "tests/nf-test-out/flashlfq/QuantifiedProteins.tsv",
    "intersected_searches": f"{pth}/Combined/unified_groups.tsv",
    "unmatched_peptides": f"{pth}/Quantify/Unmatched/unmatched_peptides.tsv",
    "open_searches": f"{pth}/Open_search/grouped_open_searches.tsv",
    "pep_threshold": 1,
    "fdr_threshold": 0.05,
}

m = mq.main(args)
