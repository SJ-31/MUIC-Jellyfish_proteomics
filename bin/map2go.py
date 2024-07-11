#!/usr/bin/env python
#
import re
import pandas as pd


def add_gos(receive, add):
    if pd.isna(receive) and pd.isna(add):
        return None
    if pd.isna(receive) and pd.notna(add):
        return add
    if pd.notna(receive) and pd.isna(add):
        return receive
    r = set(re.split("[;,]", receive))
    return ";".join(r | set(add.split(";")))


def value_from_col(df, query_column, target_column, query):
    """
    Convenience function for retrieving the corresponding entry in
    `target column` where `query_column == query`
    """
    if (try_find := df.query(f"{query_column} == '{query}'")).empty:
        return None
    return try_find[target_column].iloc[0]


def gos_from_acc(
    map_df: pd.DataFrame,
    accession: str,
    is_pfam: bool,
    interpro2go=None,
    pfam_db=None,
):
    """
    Retrieves the GO accessions corresponding to a single pfam domain
    EC accession or KEGG accession returning NA if the
    mapping fails
    """
    if pd.isna(accession):
        return accession
    gos = set()
    if is_pfam:
        # For pfam, `accession` could be received as a name that needs to
        # first be mapped to an accession.
        if find_accession := re.findall("^PF[0-9]+", accession):
            # When the accession is a proper PFAM acession code
            accession = find_accession[0].strip()
            gos.add(value_from_col(map_df, "accession", "GO", accession))
        else:
            # When the accession is just the name of the pfam
            gos.add(value_from_col(map_df, "name", "GO", accession))
            if pfam_accession := value_from_col(
                pfam_db, "name", "accession", accession
            ):
                gos.add(value_from_col(map_df, "accession", "GO", pfam_accession))
            # Try to map a pfam accession to a interpro accession, which we
            # then map to GO
            if interpro_accession := value_from_col(
                pfam_db, "name", "interpro_accession", accession
            ):
                gos.add(
                    value_from_col(interpro2go, "accession", "GO", interpro_accession)
                )
    else:
        gos.add(value_from_col(map_df, "accession", "GO", accession))
    gos = gos - {None, pd.NA}
    if gos:
        return ";".join(gos)
    return pd.NA


def go_from_row(row, pfam_map, interpro_map, ec_map, kegg_map, pfam_db):
    """
    From the corresponding mapping files, retrieve the
    GO ids of the corresponding row's PFAMs, EC and KEGG Reaction accessions
    """
    gos = set()
    spattern = "[;,]"

    if not pd.isna(row.get("PFAMs")):
        splits = re.split(spattern, row["PFAMs"])
        for s in splits:
            gos.add(
                gos_from_acc(
                    is_pfam=True,
                    map_df=pfam_map,
                    interpro2go=interpro_map,
                    pfam_db=pfam_db,
                    accession=s,
                )
            )
    for db, map in zip(["EC", "KEGG_Reaction"], [ec_map, kegg_map]):
        if db in row.index:
            if not pd.isna(row.get(db)):
                splits = re.split(spattern, row[db])
                for s in splits:
                    gos.add(gos_from_acc(map, s, is_pfam=False))
    gos = gos - {None, pd.NA}
    if gos:
        return ";".join(gos)
    return None


def parse_go_mapping(path, has_name):
    """
    Parse go mapping files (interpro2go, pfam2go, kegg_reaction2go, and ec2go)
    file into dataframe
    The behavior is different for latter 2 because their identifiers don't have
    any names as well, just accessions
    """
    if has_name:
        mapping_df = {"accession": [], "name": [], "GO": []}
    else:
        mapping_df = {"accession": [], "GO": []}
    with open(path, "r") as g:
        for line in g.readlines():
            if line.startswith("!"):
                continue
            splits = [r.strip() for r in re.split("[>;]", line)]
            acc = splits[0]
            if has_name:
                acc_split = acc.split(" ")
                acc = acc_split[0]
                if len(acc_split) > 1:
                    name = acc_split[1]
                else:
                    name = "NA"
            acc = re.sub(".*:", "", acc)
            go_id = splits[-1]
            if has_name:
                mapping_df["name"].append(name)
            mapping_df["accession"].append(acc)
            mapping_df["GO"].append(go_id)
    if has_name:
        group_cols = ["accession", "name"]
    else:
        group_cols = "accession"
    return (
        pd.DataFrame(mapping_df)
        .groupby(group_cols)
        .apply(lambda x: ";".join(x["GO"].to_list()))
        .to_frame()
        .reset_index()
        .rename({0: "GO"}, axis="columns")
    )


def map_all_db(
    pfam_db_path, p2g_path, i2g_path, ec2g_path, k2g_path, to_annotate: pd.DataFrame
):
    """
    Map entries in the pfam database to go accession numbers, using their names
    """
    print(f"Reading pfam entries from path {pfam_db_path}")
    print(f"Reading pfam2go file from path {p2g_path}")
    print(f"Reading ec2go file from path {ec2g_path}")
    print(f"Reading kegg_reaction2go file from path {k2g_path}")
    print(f"Reading interpro2go file from path {i2g_path}")
    if "GO" not in to_annotate.columns:
        return to_annotate
    pfam_db = pd.read_csv(pfam_db_path, sep="\t")
    interpro2go = parse_go_mapping(i2g_path, has_name=True)
    pfam2go = parse_go_mapping(p2g_path, has_name=True)
    ec2go = parse_go_mapping(ec2g_path, has_name=False)
    kegg_reaction2go = parse_go_mapping(k2g_path, has_name=False)
    retrieved_gos = to_annotate.apply(
        go_from_row,
        pfam_map=pfam2go,
        interpro_map=interpro2go,
        pfam_db=pfam_db,
        kegg_map=kegg_reaction2go,
        ec_map=ec2go,
        axis=1,
    )
    join_up = to_annotate["GO"].combine(retrieved_gos, add_gos)
    to_annotate["GO"] = join_up
    return to_annotate
