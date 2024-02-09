#!/usr/bin/env python
"""
Script for downloading proteins of other taxa used in comparison
"""
import re
import pandas as pd
import functools
from urllib.error import HTTPError
from io import StringIO
import requests
from requests.adapters import HTTPAdapter, Retry

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
DATABASE = "uniprotkb"
FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "length",
    "ec",
    "go",
    "go_id",
    "xref_panther",
    "xref_interpro",
    "xref_pfam",
    "xref_eggnog",
    "reviewed",
    "lineage",
]
FIELD_STR = functools.reduce(lambda x, y: f"{x}%2C{y}", FIELDS)
COMPRESSION = "false"
FORMAT = "tsv"
SIZE = 500
re_next_link = re.compile(r'<(.+)>; rel="next"')

retries = Retry(
    total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504]
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def getNextLink(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def checkResponse(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise HTTPError("API request failed!")


def getBatch(batch_url):
    """
    Generator
    """
    while batch_url:
        response = session.get(batch_url)
        checkResponse(response)
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = getNextLink(response.headers)


def tsvToFrame(response):
    return pd.read_csv(
        StringIO("\n".join(response.text.splitlines())), sep="\t"
    )


def processQuery(query) -> pd.DataFrame:
    all_results = pd.DataFrame()
    url_string = (
        f"{API_URL}/{DATABASE}/search?compressed={COMPRESSION}&fields="
        f"{FIELD_STR}&format={FORMAT}&query={query}&size={SIZE}"
    )
    for batch, total in getBatch(url_string):
        all_results = pd.concat([all_results, tsvToFrame(batch)])
    return all_results


def constructQuery(keywords, taxon_id) -> str:
    or_words = "+OR+".join(keywords)
    return f"%28%28{or_words}%29+AND+%28taxonomy_id%3A{taxon_id}%29%29"


def main(output_folder):
    taxa = {
        "Chelicerata": 6843,
        "Myriapoda": 61985,
        "Hexapoda": 6960,
        "Mammalia": 40674,
        "Serpentes": 8570,
        "Mollusca": 6447,
        "Actinopterygii": 7898,
    }

    venom_keywords = [
        "Hemolysin",
        "Venom",
        "Toxin",
        "Poison",
        "Metalloprotease",
        "Phospholipase",
    ]
    taxa_dfs = {}
    for taxon, taxon_id in taxa.items():
        search = constructQuery(venom_keywords, taxon_id)
        cur_dict = {"all": processQuery(search)}
        cur_dict["all"]["GO_IDs"] = cur_dict["all"]["Gene Ontology IDs"].apply(
            lambda x: x if pd.isna(x) else x.replace(" ", "")
        )
        cur_dict["reviewed"] = cur_dict["all"].query("Reviewed == 'reviewed'")
        taxa_dfs[taxon] = cur_dict
        cur_dict["all"].to_csv(
            f"{output_folder}/{taxon}_all.tsv",
            sep="\t",
            na_rep="NA",
            index=False,
        )
        # Filters by reviewed entries only
        cur_dict["reviewed"].to_csv(
            f"{output_folder}/{taxon}_reviewed.tsv",
            sep="\t",
            na_rep="NA",
            index=False,
        )
    return taxa_dfs


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_folder", required=True)
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    main(args["output_folder"])
