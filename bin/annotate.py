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
#

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


def parse_db(name: str, db_list: list):
    # Collect ids from a list of databases
    anno: list = []
    evi: list = []
    for db in db_list:
        db_id = db["id"]
        cur_property = db["properties"][0]["value"]
        if name == "GO":
            evidence = db["properties"][1]["value"]
            evi.append(evidence)
        if name in {"GO", "PANTHER", "Pfam"}:
            anno.append(f"{db_id}_{cur_property}")
        elif name == "InterPro":
            anno.append(db_id)
            evi.append(cur_property)
        elif name in {"KEGG", "OrthoDb", "InterPro"}:
            anno.append(db_id)
    if name == "InterPro":
        return {name: ",".join(anno), "description": ",".join(evi)}
    if name == "GO":
        return {name: ",".join(anno), "evidence": ",".join(evi)}
    return {name: ",".join(anno), "evidence": None}


def from_db(name, database_list):
    """Combine database_search and parse_db function."""
    if db := database_search(name, database_list):
        return parse_db(name, db)
    if name == "GO":
        return {"GO": None, "evidence": None}
    if name == "InterPro":
        return {"InterPro": None, "description": None}
    return {name: None, "evidence": None}


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
        "GO_evidence": [],
        "KEGG_Genes": [],
        "PFAMs": [],
        "PANTHER": [],
        "interpro_accession": [],
        "interpro_description": [],
        "length": [],
        "molWeight": [],
        "sequence": [],
    }
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results_dict: dict = get_id_mapping_results_search(link, "json")
    else:
        return pd.DataFrame(), []
    for result in results_dict["results"]:
        anno_dict["query"].append(result["from"])
        current = result["to"]
        anno_dict["UniProtKB_ID"].append(current["uniProtkbId"])
        if organism := current.get("organism"):
            anno_dict["organism"].append(organism["scientificName"])
            anno_dict["lineage"].append(",".join(organism["lineage"]))
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
            interpro = from_db("InterPro", databases)
            anno_dict["interpro_accession"].append(interpro["InterPro"])
            anno_dict["interpro_description"].append(interpro["description"])
            anno_dict["KEGG_Genes"].append(from_db("KEGG", databases)["KEGG"])
            # The 3 letters prefixed to kegg entries here are
            #   the kegg organism codes
            anno_dict["PFAMs"].append(from_db("Pfam", databases)["Pfam"])
            go = from_db("GO", databases)
            anno_dict["GO"].append(go["GO"])
            anno_dict["GO_evidence"].append(go["evidence"])
            anno_dict["PANTHER"].append(
                from_db("PANTHER", databases)["PANTHER"]
            )
        else:
            for key in [
                "KEGG_Genes",
                "GO",
                "interpro_accession",
                "interpro_description",
                "PFAMs",
                "PANTHER",
                "GO_evidence",
            ]:
                anno_dict[key].append(None)
    if failedIds := results_dict.get("failedIds"):
        failed = [f for f in failedIds]
        return anno_dict, failed
    return anno_dict, []


def idFromHeader(row):
    header = row["header"]
    cur_id = row["ProteinId"]
    if "|" in header and (find := re.search(r"\|(.*)\|", header)):
        return [cur_id, find.groups()[0]]
    elif "." in header:
        return [cur_id, header.split(" ")[0]]
    return cur_id, "NONE"


def writeFasta(needs_annotating: pd.DataFrame, file_name: str):
    query_string = ""
    for row in needs_annotating.iterrows():
        query_string = (
            query_string + f'>{row[1]["ProteinId"]}\n{row[1]["seq"]}\n'
        )
    with open(file_name, "w") as a:
        a.write(query_string)


def anno(args: dict):
    to_map = pd.read_csv(args["input"], sep="\t")
    dbIds_ids = to_map.apply(idFromHeader, axis=1).dropna()
    all_ids_mapping = pd.DataFrame(
        {
            "dbId": dbIds_ids.apply(lambda x: x[1]),
            "ProteinId": dbIds_ids.apply(lambda x: x[0]),
        }
    )
    ids = all_ids_mapping["dbId"]

    # Try to map proteins in NCBI to UniProt proteins. Multiple databases
    # required.
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
    # With every pass through a database, the list (set) of ncbi_ids will be
    # reduced as the ids are mapped by UniProt
    while len(ncbi_ids) > 0 and db_counter < 3:
        ncbi_query = map_list(list(ncbi_ids), databases[db_counter])
        anno_df = pd.DataFrame(ncbi_query[0])
        ncbi_mapped = pd.concat([ncbi_mapped, anno_df])
        ncbi_ids = ncbi_ids - set(anno_df["query"])
        db_counter += 1
    needs_annotating = pd.DataFrame()
    needs_annotating_ids = pd.Series()

    uniprot_entries = pd.DataFrame(map_list(uniprot_ids, "UniProtKB")[0])

    # (Optional) Get entries that were mapped and found in the uniprot and
    # ncbi databases, but have no annotation information
    if args["annotate_extra"]:
        no_annotation_criteria = "GO.isna() & PANTHER.isna()"
        no_annotation_criteria_inv = "GO.notna() | PANTHER.notna()"
        no_info = pd.concat(
            [
                uniprot_entries.query(no_annotation_criteria),
                ncbi_mapped.query(no_annotation_criteria),
            ]
        )["query"]
        unannotated_set = set(no_info) | ncbi_ids
        uniprot_entries = uniprot_entries.query(no_annotation_criteria_inv)
        ncbi_mapped = ncbi_mapped.query(no_annotation_criteria_inv)
    else:
        unannotated_set = ncbi_ids

    # Get any remaining unannotated sequences, ids that couldn't be mapped
    if len(ncbi_ids) > 0:
        needs_annotating = all_ids_mapping[
            all_ids_mapping["dbId"].isin(pd.Series(list(unannotated_set)))
        ]
        needs_annotating = to_map.merge(needs_annotating, on="ProteinId").drop(
            "dbId", axis="columns"
        )
        needs_annotating_headers = needs_annotating.apply(
            idFromHeader, axis=1
        ).apply(lambda x: x[1])
        needs_annotating["NCBI_ID"] = needs_annotating_headers
        needs_annotating_ids = needs_annotating["ProteinId"]
        writeFasta(needs_annotating, "needs_annotating.fasta")
        # Any remaining unannotated sequences will be extracted and sent to
        # eggnog mapper for annotation

    final = (
        pd.concat([ncbi_mapped, uniprot_entries])
        .merge(all_ids_mapping, left_on="query", right_on="dbId")
        .merge(to_map, on="ProteinId")
        .drop(["dbId", "length_y", "molWeight", "sequence"], axis="columns")
        .rename({"query": "NCBI_ID", "length_x": "length"}, axis="columns")
    )
    final = pd.concat([final, needs_annotating])
    if not needs_annotating_ids.empty:
        final[final["ProteinId"].isin(needs_annotating_ids)].drop(
            "PFAMs", axis="columns"
        ).to_csv("needs_annotating.tsv", sep="\t", index=False, na_rep="NA")
    return final


def reformat_float(x):
    if x == "-":
        return None
    return float(x)


def mergeAnnotatedEggnog(args):
    print("-" * 10)
    print("MERGING EGGNOG")
    print("-" * 10)
    eggnog_anno = pd.read_csv(args["eggnog_tsv"], sep="\t").drop(
        ["organism", "lineage", "GO_evidence", "KEGG_Genes", "PANTHER"],
        axis="columns",
    )
    eggnog_anno["NCBI_ID"] = eggnog_anno.apply(idFromHeader, axis=1).apply(
        lambda x: x[1]
    )
    anno = pd.read_csv(args["more_anno"], sep="\t")
    unwanted = set(
        eggnog_anno.columns[eggnog_anno.columns.isin(anno.columns)]
    ) - {"ProteinId"}
    get_eggnog = anno.drop(list(unwanted), axis="columns").merge(
        eggnog_anno, on="ProteinId"
    )
    get_eggnog["Anno_method"] = "eggNOG"
    already_matched = anno[~anno["ProteinId"].isin(get_eggnog["ProteinId"])]
    merged = pd.concat([already_matched, get_eggnog])
    print("-" * 10)
    print("MERGING EGGNOG DONE")
    print("-" * 10)
    return merged


def mergeAnnotatedInterpro(args):
    print("-" * 10)
    print("MERGING INTERPRO")
    print("-" * 10)
    command = "Rscript -e 'source(\"{r_source}/sort_interpro.r\")'".format(
        r_source=args["r_source"]
    )
    command = command + " -e 'df <- read_tsv(\"{input}\")'".format(
        input=args["input"]
    )
    command = command + " -e 'cleaned <- clean_annotations(df)'"
    command = command + " -e 'write_tsv(cleaned, \"sorted.tsv\")'"
    subprocess.run(command, shell=True)
    ip_sorted = pd.read_csv("./sorted.tsv", sep="\t")
    anno = pd.read_csv(args["more_anno"], sep="\t", low_memory=False)
    with open(args["interpro_query"], "r") as u:
        query_headers = [
            u.strip().replace(">", "") for u in u.readlines() if ">" in u
        ]
        query_headers = pd.Series(query_headers)
    unannotated = anno[anno["ProteinId"].isin(query_headers)]
    # Merge with metadata and drop redundant columns
    joined = ip_sorted.merge(
        unannotated, left_on="query", right_on="ProteinId"
    )
    duplicate_cols_x = joined.columns[joined.columns.str.contains("_x")]
    underscore_removed = [re.sub("_[xy]", "", c) for c in duplicate_cols_x]
    joined = joined.rename(
        dict(zip(duplicate_cols_x, underscore_removed)), axis="columns"
    ).drop([re.sub("_x", "_y", c) for c in duplicate_cols_x], axis="columns")

    annotated = anno[~(anno["ProteinId"].isin(query_headers))]
    still_left = unannotated[
        ~unannotated["ProteinId"].isin(joined["ProteinId"])
    ]
    joined["Anno_method"] = "interpro"
    # Remaining proteins that are still unannotated
    writeFasta(still_left, "still_unannotated.fasta")
    final = pd.concat([joined, annotated])
    if unwanted := set(final.columns) & {"header_x", "header_y"}:
        for u in unwanted:
            final.drop(u, axis="columns", inplace=True)
    print("-" * 10)
    print("MERGING INTERPRO DONE")
    print("-" * 10)
    return final


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--merge_interpro", action="store_true")
    parser.add_argument("-a", "--annotate_extra", action="store_true")
    parser.add_argument("-e", "--merge_eggnog", action="store_true")
    parser.add_argument("-o", "--output")
    parser.add_argument("-q", "--interpro_query")
    parser.add_argument("-r", "--r_source")
    parser.add_argument("-i", "--input")
    parser.add_argument("-m", "--more_anno")
    parser.add_argument("--eggnog_tsv")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    if args["merge_eggnog"]:
        m = mergeAnnotatedEggnog(args)
    elif args["merge_interpro"]:
        m = mergeAnnotatedInterpro(args)
    else:
        m = anno(args)
    m.to_csv(args["output"], sep="\t", index=False, na_rep="NA")
