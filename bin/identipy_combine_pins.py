#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
from pathlib import Path


def get_row(row_line, prot_col, path, excess):
    fields = row_line
    fields[0] = f"{path}.{fields[0]}"
    if excess:
        return fields[:23] + [fields[25]] + ['\t'.join(fields[prot_col:])]
    return fields[:prot_col] + ['\t'.join(fields[prot_col:])]


pins = [p for p in Path(".").glob("*.pin")]


def read_pin(path):
    name = path.stem
    with open(path, "r") as r:
        content = [line.split("\t") for line in r.read().split("\n")][:-1]
    header = content[0]
    if len(header) > 25:
        header = header[:23] + header[25:]
        lines = [get_row(l, 26, name, True) for l in content[1:]]
    else:
        lines = [get_row(l, 24, name, False) for l in content[1:]]
    df = pd.DataFrame(data=np.array(lines), columns=header)
    return (df)


all_pins = [read_pin(p) for p in pins]
all_pin_df = pd.concat(all_pins, axis=0)
all_pin_df.to_csv(sys.argv[1], sep="\t", index=False)
