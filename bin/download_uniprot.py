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
COMPRESSION = "false"
FORMAT = "tsv"
SIZE = 500
re_next_link = re.compile(r'<(.+)>; rel="next"')

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def getFieldStr(fields: list):
    return functools.reduce(lambda x, y: f"{x}%2C{y}", fields)


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
    return pd.read_csv(StringIO("\n".join(response.text.splitlines())), sep="\t")


def processQueryFasta(query: str) -> str:
    all_results = []
    url_string = (
        f"{API_URL}/{DATABASE}/search?compressed={COMPRESSION}"
        f"&format=fasta&query={query}&size={SIZE}"
    )
    for batch, total in getBatch(url_string):
        all_results.extend(batch.text.splitlines())
    return "\n".join(all_results)


def processQueryTSV(query: str, fields: str) -> pd.DataFrame:
    all_results = pd.DataFrame()
    url_string = (
        f"{API_URL}/{DATABASE}/search?compressed={COMPRESSION}&fields="
        f"{fields}&format=tsv&query={query}&size={SIZE}"
    )
    for batch, total in getBatch(url_string):
        all_results = pd.concat([all_results, tsvToFrame(batch)])
    return all_results


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fields", action="extend", nargs="+", type=str)
    parser.add_argument("--format", type=str)
    parser.add_argument("--query", type=str)
    parser.add_argument("--output", type=str)
    args = vars(parser.parse_args())  # convert to dict
    return args


# List of fields can be found here
# https://www.uniprot.org/help/return_fields

if __name__ == "__main__":
    args = parse_args()
    if args["format"] == "tsv":
        field_str = getFieldStr(args["fields"])
        tsv = processQueryTSV(args["query"], field_str)
        tsv.to_csv(args["output"], sep="\t", index=False)
    if args["format"] == "fasta":
        fasta = processQueryFasta(args["query"])
        with open(args["output"], "w") as w:
            w.write(fasta)
