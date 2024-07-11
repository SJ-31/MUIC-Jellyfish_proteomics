#!/usr/bin/env python
"""
Download matches and other features from InterPro for a given UniProt accession

Requires python >= 3.6

Example of running command:
$ python fetch-protein-matches.py UNIPROT-ACCESSION
"""

import sys
import json
import ssl
import polars as pl
from urllib import request
from urllib.error import HTTPError
from time import sleep
from urllib.request import urlopen


def from_uniprot(query):
    query = "protein/UniProt/P50876"
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/{query}"
    with urlopen(url) as res:
        data = json.loads(res.read().decode("utf-8"))
    return data


def json2df(payload) -> pl.DataFrame:
    """
    Helper function for parsing pfam json response into a dataframe
    """
    pfam_dict = {
        "accession": [],
        "name": [],
        "type": [],
        "interpro_accession": [],
    }
    for result in payload["results"]:
        print(result)
        md = result["metadata"]
        pfam_dict["accession"].append(md["accession"])
        pfam_dict["name"].append(md["name"])
        pfam_dict["type"].append(md["type"])
        pfam_dict["interpro_accession"].append(md["integrated"])
    pfam_df = pl.DataFrame(pfam_dict)
    return pfam_df


def get_all_entries(url, json_parse_fn):
    """
    Generic function for retrieving entry data
    """
    next_page = url
    context = ssl._create_unverified_context()
    dfs = []
    attempts = 0
    while next_page:
        try:
            req = request.Request(next_page, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)
            # If the API times out due a long running query
            if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
            elif res.status == 204:
                # no data so leave loop
                break
            payload = json.loads(res.read().decode())
            next_page = payload["next"]
            attempts = 0
            dfs.append(json_parse_fn(payload))
        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                # If there is a different HTTP error, it wil re-try 3 times
                # before failing
                if attempts < 3:
                    attempts += 1
                    sleep(61)
                    continue
                else:
                    sys.stderr.write("LAST URL: " + next_page)
                    raise e
    return pl.concat(dfs)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--get_pfam")
    parser.add_argument("--get_panther")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    arguments = parse_args()
    # will conflict with reticulate otherwise
    if arguments["get_pfam"]:
        pfam_url = "https://www.ebi.ac.uk:443/interpro/api/entry/pfam/?page_size=200"
        all_pfam: pl.DataFrame = get_all_entries(pfam_url, json2df)
        all_pfam.write_csv(arguments["get_pfam"], separator="\t")
    elif arguments["get_panther"]:
        panther_url = (
            "https://www.ebi.ac.uk:443/interpro/api/entry/panther/?page_size=200"
        )
        all_pfam: pl.DataFrame = get_all_entries(panther_url, json2df)
        all_pfam.write_csv(arguments["get_panther"], separator="\t")
