import polars as pl
import sys
import os

os.chdir("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/docs")
sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import asyncio
import download_uniprot as du


def get_identifiers(df, pool_name):
    df = df.filter(pl.col("Pool name") == pool_name)
    return list(
        {
            u
            for u_str in df.filter(pl.col("Unique identifiers").is_not_null())[
                "Unique identifiers"
            ]
            for u in u_str.split(";")
        }
    )


sample = pl.read_csv(
    "../results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv",
    separator="\t",
    null_values="NA",
)

MEAN_PEPTIDE_LENGTH = (
    pl.Series(
        list(
            {
                p
                for peptide_str in sample.filter(
                    pl.col("unique_peptides").is_not_null()
                )["unique_peptides"]
                for p in peptide_str.split(";")
            }
        )
    )
    .str.len_chars()
    .mean()
)
print(f"Mean peptide length in sample: {MEAN_PEPTIDE_LENGTH}")


def get_pool_stats(df: pl.DataFrame) -> pl.DataFrame:
    return df.group_by("Pool name").agg(
        pl.col("Sequence").str.len_chars().mean().alias("Mean sequence length"),
        pl.col("Sequence").str.len_chars().median().alias("Median sequence length"),
        pl.col("Gene names").null_count().alias("N missing gene names"),
        pl.col("Unique identifiers").null_count().alias("N missing identifiers"),
    )


def choose_pool(df: pl.DataFrame) -> str:
    def print_info(df2):
        print(f"Best choice: {df2['Pool name'].item()}")
        print(f"\tN missing identifiers = {df2['N missing identifiers'].item()}")
        print(f"\tN missing gene names = {df2['N missing gene names'].item()}")
        print(f"\tMean sequence length = {df2['Mean sequence length'].item()}")
        return df2["Pool name"].item()

    stats = get_pool_stats(df).sort("N missing identifiers")
    choices = stats.filter(
        pl.col("N missing identifiers") == min(stats["N missing identifiers"])
    )
    if choices.shape[0] == 1:
        print("Found best by identifiers")
        return print_info(choices)
    choices = choices.with_columns(
        distance_from_sample_mean=(
            pl.col("Mean sequence length") - MEAN_PEPTIDE_LENGTH
        ).abs()
    )
    choices = choices.filter(
        pl.col("distance_from_sample_mean") == min(choices["distance_from_sample_mean"])
    )
    if choices.shape[0] == 1:
        print("Found best by distance from mean peptide length")
        return print_info(choices)
    choices = choices.sample(n=1)
    print(
        "Cannot decide by mean peptide length or missing identifiers, taking random..."
    )
    return print_info(choices)


missing_gene_set = pl.read_csv(
    "../data/benchmark/ProteomeTools/MissingGeneSet.tsv", separator="\t"
)
srmatlas = pl.read_csv(
    "../data/benchmark/ProteomeTools/SRMAtlas Set_.tsv", separator="\t"
)
proteotypic = pl.read_csv(
    "../data/benchmark/ProteomeTools/Proteotypic Set_.tsv", separator="\t"
)

print("\nProteotypic")
pr = choose_pool(proteotypic)
print("\nMissing gene")
mg = choose_pool(missing_gene_set)
print("\nSRMAtlas")
sa = choose_pool(srmatlas)

all_ids = list(
    set(
        get_identifiers(proteotypic, pr)
        + get_identifiers(missing_gene_set, mg)
        + get_identifiers(srmatlas, sa)
    )
)


def make_id_query(id) -> str:
    return f"%28%28accession%3A{id}%29%29"  # Technically accession, not id



results = {}
for id in all_ids:
    results[id] = du.process_query_fasta(make_id_query(id))

def prepend_to_fasta(fasta_str, acc):
    return f">ACC:{acc} {fasta_str.replace(">", "")}"

failed_ids = list(filter(lambda x: not results[x], results))
# Some accessions failed because they have been merged into other entries

ids: list = list(results.keys())
for i in ids:
    if i in failed_ids:
        del results[i]
        continue
    results[i] = prepend_to_fasta(results[i], i)


with open(f"../data/benchmark/ProteomeTools/failed_ids.txt", "w") as w:
    w.write("\n".join(failed_ids))

filename = f"{mg}_{sa}_{pr}.fasta"
with open(f"../data/benchmark/ProteomeTools/{filename}", "w") as w:
    w.write("\n".join(list(results.values())))
