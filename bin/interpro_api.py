#!/usr/bin/env ipython
"""
Download matches and other features from InterPro for a given UniProt accession

Requires python >= 3.6

Example of running command:
$ python fetch-protein-matches.py UNIPROT-ACCESSION
"""

import json
import sys
from urllib.error import HTTPError
from urllib.request import urlopen

def from_uniprot(query):
    query = "protein/UniProt/P50876"
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/{query}"
    with urlopen(url) as res:
        data = json.loads(res.read().decode("utf-8"))
    return data


dbs = ["interpro", "smart", "pfam"]
query = "entry/interpro/IPR011106"
api_url = "https://www.ebi.ac.uk/interpro/api"
url = f"{api_url}/{query}"
with urlopen(url) as res:
    data = json.loads(res.read().decode("utf-8"))
