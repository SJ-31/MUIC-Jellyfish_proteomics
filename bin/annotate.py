#!/usr/bin/env python
import subprocess
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

retries = Retry(
    total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504]
)
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
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
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
            return bool(j.get("results") or j.get("failedIds"))


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
                line
                for line in decompressed.decode("utf-8").split("\n")
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
    compressed = (
        query["compressed"][0].lower() == "true"
        if "compressed" in query
        else False
    )
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
    if db := database_search(name, database_list):
        return parse_db(name, db)
    return "None"


def map_list(id_list, origin_db: str):
    job_id = submit_id_mapping(
        from_db=origin_db, to_db="UniProtKB", ids=id_list
    )
    anno_dict: dict = {
        "query": [],
        "UniProtKB_ID": [],
        "organism": [],
        "lineage": [],
        "GO": [],
        "KEGG_Genes": [],
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
        return (pd.DataFrame(), [])
    for result in results_dict["results"]:
        anno_dict["query"].append(result["from"])
        current = result["to"]
        anno_dict["UniProtKB_ID"].append(current["uniProtkbId"])
        if organism := current.get("organism"):
            anno_dict["organism"].append(organism["scientificName"])
            anno_dict["lineage"].append(";".join(organism["lineage"]))
        else:
            anno_dict["organism"].append("unknown")
            anno_dict["lineage"].append("unknown")
        if sequence := current.get("sequence"):
            anno_dict["length"].append(sequence["length"])
            anno_dict["molWeight"].append(sequence["molWeight"])
            anno_dict["sequence"].append(sequence["value"])
        else:
            for key in ["length", "molWeight", "sequence"]:
                anno_dict[key].append("unknown")
        if databases := current.get("uniProtKBCrossReferences"):
            anno_dict["KEGG_Genes"].append(from_db("KEGG", databases))
            # The 3 letters prefixed to kegg entries here are
            #   the kegg organism codes
            anno_dict["GO"].append(from_db("GO", databases))
            anno_dict["PANTHER"].append(from_db("PANTHER", databases))
            anno_dict["OrthoDb"].append(from_db("OrthoDb", databases))
        else:
            for key in ["KEGG_Genes", "GO", "PANTHER", "OrthoDb"]:
                anno_dict[key].append("NA")
    if failedIds := results_dict.get("failedIds"):
        failed = [f for f in failedIds]
        return (anno_dict, failed)
    return (anno_dict, [])


def id_from_header(row):
    header = row["header"]
    id = row["ProteinId"]
    if "|" in header and (find := re.search(r"\|(.*)\|", header)):
        return [id, find.groups()[0]]
    elif "." in header:
        return [id, header.split(" ")[0]]
    return (id, "NONE")


def write_fasta(needs_annotating: pd.DataFrame, file_name: str):
    query_string = ""
    for row in needs_annotating.iterrows():
        query_string = (
            query_string + f'>{row[1]["ProteinId"]}\n{row[1]["seq"]}\n'
        )
    with open(file_name, "w") as a:
        a.write(query_string)


def main(args: dict):
    to_map = pd.read_csv(args["input"], sep="\t")
    dbIds_ids = to_map.apply(id_from_header, axis=1).dropna()
    id_map = pd.DataFrame(
        {
            "dbId": dbIds_ids.apply(lambda x: x[1]),
            "ProteinId": dbIds_ids.apply(lambda x: x[0]),
        }
    )
    ids = id_map["dbId"]

    # Try to map proteins in NCBI to UniProt proteins. Multiple databases
    # required
    ncbi_ids = ids.where(ids.str.contains("\\."))
    ncbi_ids = set(ncbi_ids.dropna())
    uniprot_ids = ids[~(ids.isin(ncbi_ids))]
    ncbi_mapped = pd.DataFrame()
    databases = [
        "RefSeq_Protein",
        "EMBL-GenBank-DDBJ",
        "EMBL-GenBank-DDBJ_CDS",
    ]
    db_counter = 0
    while len(ncbi_ids) > 0 and db_counter < 3:
        ncbi_query = map_list(list(ncbi_ids), databases[db_counter])
        anno_df = pd.DataFrame(ncbi_query[0])
        ncbi_mapped = pd.concat([ncbi_mapped, anno_df])
        ncbi_ids = ncbi_ids - set(anno_df["query"])
        db_counter += 1
    needs_annotating = pd.DataFrame()
    needs_annotating_headers = pd.Series()

    # Get any remaining unannotated sequences
    if len(ncbi_ids) > 0:
        needs_annotating = id_map[
            id_map["dbId"].isin(pd.Series(list(ncbi_ids)))
        ]
        needs_annotating = to_map.merge(needs_annotating, on="ProteinId")
        needs_annotating_headers = needs_annotating["ProteinId"]
        write_fasta(needs_annotating, "needs_annotating.fasta")
        # Any remaining unannotated sequences will be extracted and sent to
        # interproscan for annotation
    uniprot_query = map_list(uniprot_ids, "UniProtKB")
    final = (
        pd.concat([ncbi_mapped, pd.DataFrame(uniprot_query[0])])
        .merge(id_map, left_on="query", right_on="dbId")
        .merge(to_map, on="ProteinId")
        .drop(["dbId", "length_y", "molWeight", "sequence"], axis="columns")
        .rename({"query": "NCBI_ID", "length_x": "length"}, axis="columns")
    )
    final = pd.concat([final, needs_annotating])
    anno_cols = [
        "NCBI_ID",
        "UniProtKB_ID",
        "organism",
        "lineage",
        "GO",
        "KEGG_Genes",
        "OrthoDb",
        "PANTHER",
    ]
    not_anno = list(set(final.columns) - set(anno_cols))
    if not needs_annotating_headers.empty:
        final[final["ProteinId"].isin(needs_annotating_headers)].to_csv(
            "needs_annotating.tsv", sep="\t", index=False, na_rep="NaN"
        )
    anno = final[["ProteinId", "header"] + anno_cols]
    meta = final[["ProteinId"] + not_anno]
    anno.to_csv(args["anno_tsv"], sep="\t", index=False)
    meta.to_csv(args["meta_tsv"], sep="\t", index=False)


def reformat_float(x):
    if x == "-":
        return None
    return float(x)


def merge_annotated_eggnog(args):
    eggnog_anno = pd.read_csv(args["eggnog_anno_tsv"], sep="\t")
    anno = pd.read_csv(args["anno_tsv"], sep="\t")
    get_eggnog = anno.drop("GO", axis="columns").merge(
        eggnog_anno, on="ProteinId"
    )
    already_matched = anno[~anno["ProteinId"].isin(get_eggnog["ProteinId"])]
    merged = pd.concat([already_matched, get_eggnog])
    merged.to_csv(args["anno_tsv"], sep="\t", index=False)


def merge_annotated_interpro(args):
    command = "Rscript -e 'source(\"{r_source}/sort_interpro.r\")'".format(
        r_source=args["r_source"]
    )
    command = command + " -e 'df <- read_tsv(\"{input}\")'".format(
        input=args["input"]
    )
    command = command + " -e 'cleaned <- clean_annotations(df)'"
    command = command + " -e 'write_tsv(cleaned, \"sorted.tsv\")'"
    subprocess.run(command, shell=True)
    sorted = pd.read_csv("./sorted.tsv", sep="\t")
    anno = pd.read_csv(args["anno_tsv"], sep="\t")
    with open(args["interpro_query"], "r") as u:
        query_headers = [
            u.strip().replace(">", "") for u in u.readlines() if ">" in u
        ]
        query_headers = pd.Series(query_headers)
    unannotated = anno[anno["ProteinId"].isin(query_headers)]
    joined = (
        sorted.merge(unannotated, left_on="query", right_on="ProteinId")
        .rename({"GO_x": "GO"}, axis="columns")
        .drop(["GO_y", "query"], axis="columns")
    )
    joined["organism"] = joined["header"].apply(
        lambda x: x[x.find("[") + 1 : x.find("]")] if "[" in x else "-"
    )
    annotated = anno[~(anno["ProteinId"].isin(query_headers))]
    still_left = unannotated[
        ~unannotated["ProteinId"].isin(joined["ProteinId"])
    ]
    # Remaining proteins that are still unannotated
    write_fasta(still_left, "still_needs_annotating.fasta")
    final = pd.concat([joined, annotated])
    final.to_csv(args["anno_tsv"], sep="\t", index=False)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--merge_interpro", action="store_true")
    parser.add_argument("-e", "--merge_eggnog", action="store_true")
    parser.add_argument("-m", "--meta_tsv")
    parser.add_argument("-q", "--interpro_query")
    parser.add_argument("-a", "--anno_tsv")
    parser.add_argument("-r", "--r_source")
    parser.add_argument("-i", "--input")
    parser.add_argument("--eggnog_anno_tsv")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    if args["merge_eggnog"]:
        merge_annotated_eggnog(args)
    elif args["merge_interpro"]:
        merge_annotated_interpro(args)
    else:
        main(args)
