#!/usr/bin/env ipython
pth = "./results/jellyfish/1-First_pass"
args = {
    "direct_lfq": f"{pth}/Quantify/directlfq_prot.tsv",
    "intersected_searches": f"{pth}/Combined/unified_groups.tsv",
    "open_searches": f"{pth}/Open_search/grouped_open_searches.tsv",
    "pep_threshold": 1,
    "fdr_threshold": 0.05,
}
