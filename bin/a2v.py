#!/usr/bin/env python

import re
import csv
import numpy as np
import pandas as pd


def loadEmbeddings(path, key, selected_gos):
    all_embeddings = np.load(path, allow_pickle=True)[key].item()
    unwanted = set(all_embeddings.keys()) - selected_gos
    for key in unwanted:
        del all_embeddings[key]
    return all_embeddings


def getGos(results_file):
    all_gos = set()
    with open(results_file, "r") as r:
        rd = csv.reader(r, delimiter="\t")
        header = next(rd)
        if "GO" in header:
            n = 0
            while header[n] != "GO":
                n += 1
            go_loc = n
        else:
            raise ValueError("GO header not in the input file!")
        for row in rd:
            gos = row[go_loc]
            for g in gos.split(";"):
                if g != "NA":
                    cleaned = re.findall("GO:[0-9]+", g)[0]
                    all_gos.add(cleaned)
    return all_gos
