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


def get_field_str(fields: list):
    return functools.reduce(lambda x, y: f"{x}%2C{y}", fields)


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise HTTPError("API request failed!")


def get_batch(batch_url):
    """
    Generator
    """
    while batch_url:
        response = session.get(batch_url)
        check_response(response)
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


def tsv2frame(response):
    return pd.read_csv(StringIO("\n".join(response.text.splitlines())), sep="\t")


def process_query_fasta(query: str) -> str:
    all_results = []
    url_string = (
        f"{API_URL}/{DATABASE}/search?compressed={COMPRESSION}"
        f"&format=fasta&query={query}&size={SIZE}"
    )
    for batch, total in get_batch(url_string):
        all_results.extend(batch.text.splitlines())
    return "\n".join(all_results)


def process_query_tsv(query: str, fields: str) -> pd.DataFrame:
    all_results = pd.DataFrame()
    url_string = (
        f"{API_URL}/{DATABASE}/search?compressed={COMPRESSION}&fields="
        f"{fields}&format=tsv&query={query}&size={SIZE}"
    )
    for batch, total in get_batch(url_string):
        all_results = pd.concat([all_results, tsv2frame(batch)])
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
        field_str = get_field_str(args["fields"])
        tsv = process_query_tsv(args["query"], field_str)
        tsv.to_csv(args["output"], sep="\t", index=False)
    if args["format"] == "fasta":
        fasta = process_query_fasta(args["query"])
        with open(args["output"], "w") as w:
            w.write(fasta)
