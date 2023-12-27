#!/usr/bin/env ipython
import os
import importlib
import pandas as pd
import sys

os.chdir("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/pytest")
sys.path.append("../../bin")
import msgf2pin as m2p


def join_proteins(bad_line: list[str]):
    fits = bad_line[:32]
    the_rest = bad_line[32:]
    if not the_rest:
        return fits
    return fits[:-1] + [",".join([fits[-1]] + the_rest)]


dirname = "../../results/jellyfish/1-First_pass"
args = {
    "decoy_mzid": f"{dirname}/Engines/MSGF/decoys-CiCs3.mzid",
    "valid_mzid": f"{dirname}/Engines/MSGF/normal-CiCs3.mzid",
    "real": "../results/msgf2pin/sample_pin.tab",
    "output": "../results/msgf2pin/my-convert.pin",
}
compare = pd.read_table(
    args["real"], sep="\t", on_bad_lines=join_proteins, engine="python"
)

m2p.main(args)
check = pd.read_table(
    args["output"], sep="\t", on_bad_lines=join_proteins, engine="python"
)
