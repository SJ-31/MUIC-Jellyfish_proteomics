#!/usr/bin/env python
#
import polars as pl
import polars.selectors as cs


def write_fasta(headers: list[str], seqs: list[str], filename: str) -> None:
    text = "\n".join([f">{h}\n{s}" for h, s in zip(headers, seqs)])
    with open(filename, "w") as w:
        w.write(text)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast_results")
    parser.add_argument("-u", "--unknown_hits")
    parser.add_argument("-q", "--blast_query_map")
    parser.add_argument("-d", "--database_hits")
    parser.add_argument("-m", "--seq_mapping")
    parser.add_argument("-f", "--unmatched_fasta_output")
    parser.add_argument("-t", "--unmatched_output")
    parser.add_argument("-a", "--accepted_output")
    parser.add_argument("-o", "--matched_output")
    parser.add_argument("-i", "--identity_threshold", type=float)
    parser.add_argument("-p", "--pep_threshold", type=float)
    parser.add_argument("-e", "--evalue_threshold", type=float)
    args = vars(parser.parse_args())  # convert to dict
    return args


def main(args: dict):
    blast_df = pl.read_csv(args["blast_results"])

    qmap = pl.read_csv(args["blast_query_map"], separator="\t")
    group_df = pl.read_csv(args["database_hits"], separator="\t", null_values="NA")
    unknown_df = pl.read_csv(args["unknown_hits"], separator="\t", null_values="NA")

    # Join with the query map to obtain cleaned sequences
    filtered = (
        (
            blast_df.filter(
                (pl.col("pident") >= args["identity_threshold"])
                & (pl.col("evalue") <= args["evalue_threshold"])
            )
        )
        .join(qmap, on="queryId")
        .join(
            unknown_df.select(
                "ProteinId",
                "q.value",
                "posterior_error_prob",
                "ID_method",
                "inferred_by",
                "ProteinGroupId",
            ),
            on="ProteinId",
        )
    )

    # Write out the query sequences that were not matched successfully by BLAST
    unmatched = unknown_df.filter(
        (~pl.col("ProteinId").is_in(filtered["ProteinId"]))
    ).unique("seq")

    # Mark best blast hits and one-hit wonders
    find_best = []
    for _, data in filtered.group_by("queryId"):
        best = data["evalue"].arg_min()
        best_col = [str(int(i == best)) for i in range(data.shape[0])]
        new = data.with_columns(is_blast_best=pl.Series(best_col))
        find_best.append(new)
    found_best: pl.DataFrame = pl.concat(find_best)

    one_hits = found_best.group_by("subjectId").agg(pl.len())
    one_hit_dct = dict(zip(one_hits["subjectId"], one_hits["len"]))
    found_best = found_best.with_columns(
        is_blast_one_hit=pl.col("subjectId").map_elements(
            lambda x: str(int(one_hit_dct[x] == 1)), return_dtype=pl.String
        )
    )

    # Group up blast queries into single proteins
    to_keep = [
        "peptideIds",
        "ProteinId",
        "is_blast_best",
        "is_blast_one_hit",
        "posterior_error_prob",
        "q.value",
    ]
    join_exprs = [pl.col(c).list.join(";") for c in to_keep]
    split_unique = ["inferred_by", "ID_method", "ProteinGroupId"]
    split_unique_exprs = [pl.col(c).list.join(separator=";") for c in split_unique]
    found_best = (
        (
            found_best.group_by("subjectId")
            .agg(pl.col(to_keep), pl.col(split_unique).unique())
            .with_columns(join_exprs)
            .rename({"ProteinId": "MatchedPeptideIds", "subjectId": "ProteinId"})
        )
        .join(
            pl.read_csv(args["seq_mapping"], separator="\t"),
            left_on="ProteinId",
            right_on="id",
        )
        .with_columns(split_unique_exprs)
        .with_columns(
            pl.col(split_unique).map_elements(
                lambda x: ";".join(set(x.split(";"))), return_dtype=pl.String
            )
        )
    )

    # Join up blast identifications with the search engine results
    to_concat = ["peptideIds", "q.value", "posterior_error_prob"]
    concat_exprs = [
        pl.concat_str([c, f"{c}_right"], separator=";").alias(c) for c in to_concat
    ]
    was_in_previous: pl.DataFrame = (
        group_df.join(found_best, on="ProteinId")
        .with_columns(concat_exprs)
        .drop(cs.contains("_right"))
    )
    blast_id_only = found_best.filter(
        ~pl.col("ProteinId").is_in(was_in_previous["ProteinId"])
    ).select(was_in_previous.columns)

    null_cols = [
        pl.lit(None).alias(n)
        for n in set(was_in_previous.columns) - set(group_df.columns)
    ]
    others = (
        group_df.filter(~pl.col("ProteinId").is_in(was_in_previous["ProteinId"]))
        .with_columns(null_cols)
        .select(was_in_previous.columns)
    )

    final = pl.concat([was_in_previous, blast_id_only, others])
    return final, unmatched, found_best


if __name__ == "__main__":
    args = parse_args()
    f, u, a = main(args)
    f.write_csv(args["matched_output"], separator="\t", null_value="NA")
    u.write_csv(args["unmatched_output"], separator="\t", null_value="NA")
    write_fasta(u["ProteinId"], u["seq"], args["unmatched_fasta_output"])
    a.write_csv(args["accepted_output"], separator="\t", null_value="NA")
