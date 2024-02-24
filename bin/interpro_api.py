#!/usr/bin/env python
"""
Download matches and other features from InterPro for a given UniProt accession

Requires python >= 3.6

Example of running command:
$ python fetch-protein-matches.py UNIPROT-ACCESSION
"""

import sys
import re
import json
import ssl
import pandas as pd
from urllib import request
from urllib.error import HTTPError
from time import sleep
from urllib.request import urlopen


def valueFromCol(df, query_column, target_column, query):
    if (try_find := df.query(f"{query_column} == '{query}'")).empty:
        return None
    return try_find[target_column].iloc[0]


def from_uniprot(query):
    query = "protein/UniProt/P50876"
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/{query}"
    with urlopen(url) as res:
        data = json.loads(res.read().decode("utf-8"))
    return data


def makePfamDf(payload):
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
        md = result["metadata"]
        pfam_dict["accession"].append(md["accession"])
        pfam_dict["name"].append(md["name"])
        pfam_dict["type"].append(md["type"])
        pfam_dict["interpro_accession"].append(md["integrated"])
    pfam_df = pd.DataFrame(pfam_dict)
    return pfam_df


def getAllPfam():
    """
    Main function for retrieving all entries in the pfam database
    """
    BASE_URL = (
        "https://www.ebi.ac.uk:443/interpro/api/entry/pfam/?page_size=200"
    )
    next_page = BASE_URL
    context = ssl._create_unverified_context()
    dfs = []
    attempts = 0
    while next_page:
        try:
            req = request.Request(
                next_page, headers={"Accept": "application/json"}
            )
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
            dfs.append(makePfamDf(payload))
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
    return pd.concat(dfs)


def addGos(receive, add):
    if pd.isna(receive) and pd.isna(add):
        return None
    if pd.isna(receive) and pd.notna(add):
        return add
    if pd.notna(receive) and pd.isna(add):
        return receive
    return ";".join(set(receive.split(";")) | set(add.split(";")))


def parseGoMapping(path):
    """
    Parse go mapping (interpro2go, pfam2go) file into dataframe
    http://current.geneontology.org/ontology/external2go/pfam2go
    """
    p2g_df = {"accession": [], "name": [], "GO": []}
    with open(path, "r") as g:
        for line in g.readlines():
            if line.startswith("!"):
                continue
            go = line.split(";")[1].strip()
            pfams = line.split(">")[0].strip().split(" ")
            p2g_df["name"].append(" ".join(pfams[1:]))
            p2g_df["accession"].append(re.sub(".*:", "", pfams[0]))
            p2g_df["GO"].append(go)
    return (
        pd.DataFrame(p2g_df)
        .groupby(["accession", "name"])
        .apply(lambda x: ";".join(x["GO"].to_list()))
        .to_frame()
        .reset_index()
        .rename({0: "GO"}, axis="columns")
    )


def getOnePfam(pfam2go, interpro2go, pfam_db, domain_name):
    """
    Maps a single pfam domain to a GO accession, returning NA if the
    mapping fails
    """
    gos = set()
    gos.add(valueFromCol(pfam2go, "name", "GO", domain_name))
    if pfam_accession := valueFromCol(
        pfam_db, "name", "accession", domain_name
    ):
        gos.add(valueFromCol(pfam2go, "accession", "GO", pfam_accession))
    if interpro_accession := valueFromCol(
        pfam_db, "name", "interpro_accession", domain_name
    ):
        gos.add(
            valueFromCol(interpro2go, "accession", "GO", interpro_accession)
        )
    gos = gos - {None}
    if gos:
        return ";".join(gos)
    return None


def pfamsFromRow(row, pfam2go, interpro2go, pfam_db):
    gos = set()
    domain_name = row["PFAMs"]
    splits = domain_name.split(";")
    if len(splits) == 1:
        return getOnePfam(pfam2go, interpro2go, pfam_db, domain_name)
    for s in splits:
        gos.add(getOnePfam(pfam2go, interpro2go, pfam_db, s))
    gos = gos - {None}
    if gos:
        return ";".join(gos)
    return None


def mapPfams(pfam_db_path, p2g_path, i2g_path, to_annotate):
    """
    Map entries in the pfam database to go accession numbers, using
    """
    print(f"Reading pfam entries from path {pfam_db_path}")
    print(f"Reading pfam2go file from path {p2g_path}")
    print(f"Reading interpro2go entries from path {i2g_path}")
    has_pfams = to_annotate.query("~PFAMs.isna()")
    no_pfam = to_annotate.query("PFAMs.isna()")

    pfam_db = pd.read_csv(pfam_db_path, sep="\t")
    interpro2go = parseGoMapping(i2g_path)
    pfam2go = parseGoMapping(p2g_path)
    mapped_gos = has_pfams.apply(
        pfamsFromRow,
        pfam2go=pfam2go,
        interpro2go=interpro2go,
        pfam_db=pfam_db,
        axis=1,
    )
    join_up = has_pfams["GO"].combine(mapped_gos, addGos)
    (has_pfams)["GO"] = join_up
    return pd.concat([has_pfams, no_pfam])


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--get_pfam")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    arguments = parse_args()
    # will conflict with reticulate otherwise
    if arguments["get_pfam"]:
        all_pfam = getAllPfam()
        all_pfam.to_csv(
            arguments["get_pfam"], sep="\t", na_rep="NA", index=False
        )
