#!/usr/bin/env ipython
import pandas as pd
import re

file = pd.read_csv("test_casanovo.mztab", sep="\t", skiprows=61)
peps = file["sequence"]
print(peps[11040])
print(re.match("[^A-Z]*", peps[11040]))
