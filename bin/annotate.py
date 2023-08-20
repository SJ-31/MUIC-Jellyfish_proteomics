#!/usr/bin/env python

import sys
import re
import pandas as pd
import time
import json
import zlib
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
# Use UniProt's mapping API to obtain GO, KEGG and UniProt ids for
#   a genes in the NCBI's gene database
#   Input: a tsv file that contains protein fasta headers in a
#       column titled "header"

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

retries = Retry(total=5,
                backoff_factor=0.25,
                status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={
            "from": from_db,
            "to": to_db,
            "ids": ",".join(ids)
        },
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [
                line for line in decompressed.decode("utf-8").split("\n")
                if line
            ]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text

def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url, file_format):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (query["compressed"][0].lower() == "true"
                  if "compressed" in query else False)
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    return results


def database_search(get: str, databases: list) -> list:
    # check if a database ("get") is in the results
    found: list = []
    for db in databases:
        if ("database", get) in db.items():
            found.append(db)
    return found


def parse_db(name: str, db_list: list) -> str:
    # Collect ids from a list of databases
    anno: list = []
    for db in db_list:
        id = db["id"]
        property = db["properties"][0]["value"]
        if name in {"GO", "PANTHER", "Pfam"}:
            anno.append(f"{id}_{property}")
        elif name in {"KEGG", "OrthoDb"}:
            anno.append(id)
    return ";".join(anno)


def from_db(name, database_list) -> str:
    """Combine database_search and parse_db function."""
    if (db := database_search(name, database_list)):
        return parse_db(name, db)
    return "None"



def map_list(id_list, origin_db: str):
    job_id = submit_id_mapping(from_db=origin_db, to_db="UniProtKB", ids=id_list)
    anno_dict: dict = {
        "Gene_Name": [],
        "UniProtKB_ID": [],
        "organism": [],
        "lineage": [],
        "GO": [],
        "KEGG": [],
        "OrthoDb": [],
        "PANTHER": [],
        "length": [],
        "molWeight": [],
        "sequence": [],
    }
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results_dict: dict = get_id_mapping_results_search(link, "json")
    else:
        return None
    for result in results_dict["results"]:
        anno_dict["Gene_Name"].append(result["from"])
        current = result["to"]
        anno_dict["UniProtKB_ID"].append(current["uniProtkbId"])
        anno_dict["organism"].append(current["organism"]["scientificName"])
        anno_dict["lineage"].append(
            ';'.join(current["organism"]["lineage"]))
        anno_dict["length"].append(current["sequence"]["length"])
        anno_dict["molWeight"].append(current["sequence"]["molWeight"])
        anno_dict["sequence"].append(current["sequence"]["value"])
        databases = current["uniProtKBCrossReferences"]
        anno_dict["KEGG"].append(from_db("KEGG", databases))
        anno_dict["GO"].append(from_db("GO", databases))
        anno_dict["PANTHER"].append(from_db("PANTHER", databases))
        anno_dict["OrthoDb"].append(from_db("OrthoDb", databases))
    mapping = pd.DataFrame(anno_dict)
    return mapping

def id_from_header(row, denovo_list, transcriptome_list, seq_mapping):
    header = row["header"]
    id = row["ProteinId"]
    if "|" in header and (find := re.search(r'\|(.*)\|', header)):
        return find.groups()[0]
    elif "." in header:
        return header.split(" ")[0]
    elif "-DENOVO" in header:
        lookup = seq_mapping.query("id == @id")
        denovo_list.append(f'>{header}\n{lookup["seq"]}')
    elif "-TRANSCRIPTOME" in header:
        lookup = seq_mapping.query("id == @id")
        transcriptome_list.append(f'>{header}\n{lookup["seq"]}')
    return None

# Main
input = sys.argv[1]
output = sys.argv[2]
seq_mapfile = sys.argv[3]

# input = "testing.tsv"
# output = "sort_requests"
# seq_mapfile = "../ref/all_normal_mapping.tsv"

denovo_hits = []
transcriptome_hits = []
to_map = pd.read_csv(input, sep="\t")
seq_map = pd.read_csv(seq_mapfile, sep="\t")
ids = to_map.apply(id_from_header, denovo_list=denovo_hits,
                   transcriptome_list=transcriptome_hits,
                   seq_mapping=seq_map,
                   axis=1).dropna()
ncbi_ids = ids.where(ids.str.contains("\\.")).dropna()
# Use RefSeq_Protein
uniprot_ids = ids[~(ids.isin(ncbi_ids))]
# Use UniProtKB

final = pd.concat([map_list(ncbi_ids, "RefSeq_Protein"),
           map_list(uniprot_ids, "UniProtKB")])
final.to_csv(output, index=False)
for fasta, lst in zip(["denovo_hits", "transcriptome_hits"],
                      [denovo_hits, transcriptome_hits]):
    with open(f"{fasta}.fasta", "w", "utf-8") as f:
        f.write("".join(lst))
