#!/usr/bin/env ipython
#
import pandas as pd
import sys

pth = "./results/jellyfish/1-First_pass"
sys.path.append("./bin")
import merge_quantifications as mq
import write_quant as wq

args = {
    "dlfq": f"{pth}/Quantify/directlfq_prot.tsv",
    "flfq": "tests/nf-test-out/flashlfq/QuantifiedProteins.tsv",
    "intersected_searches": f"{pth}/Combined/unified_groups.tsv",
    "unmatched_peptides": f"{pth}/Quantify/Unmatched/unmatched_peptides.tsv",
    "open_searches": f"{pth}/Open_search/grouped_open_searches.tsv",
}

m = mq.main(args)
flfq = wq.read_flashlfq(args["flfq"])
dlfq = wq.read_directlfq(args["dlfq"])
