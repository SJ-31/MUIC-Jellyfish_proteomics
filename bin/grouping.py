import polars as pl
import sys
import polars as pl
import re

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
tf = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/toxin_groups.tsv"
d = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"
pfam_map_file = (
    "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/Pfam-A.clans.tsv"
)

REGEXES = {
    "pfam": re.compile("(PF.*?)_"),
    "uniprot": re.compile("\\S*\\|\\S*\\|\\S* (.*) OS="),
    "ncbi": re.compile("\\S* (.*) \\["),
}

ANNO_COLS = [
    "PFAMs",
    "header",
    "organism",
    "PANTHER",
    "eggNOG_description",
    "interpro_description",
    "interpro_accession",
    "eggNOG_OGs",
    "GO_IDs",
]
toxin_map = pl.read_csv(tf, separator="\t")
pfam_map = pl.read_csv(pfam_map_file, separator="\t")


def get_from_df(key, key_col, value_col, df, filter_criteria: tuple = ()):
    """
    Retrieve a value from a DataFrame based on a key and optional filter criteria.

    Args:
        key (Any): The key to search for in the DataFrame.
        key_col (str): The column name in the DataFrame where the key is located.
        value_col (str): The column name in the DataFrame from which the value should be retrieved.
        df (polars.DataFrame): The DataFrame to search within.
        filter_criteria (tuple, optional): A tuple containing a column name and value to filter the DataFrame before searching. Defaults to an empty tuple.

    Returns:
        Any: The value from the specified column corresponding to the key, or None if the key is not found.

    Raises:
        ValueError: If multiple entries in the DataFrame match the key.
    """
    if filter_criteria:
        df = df.filter(pl.col(filter_criteria[0]) == filter_criteria[1])
    try_find = df.filter(pl.col(key_col) == key)
    if try_find.is_empty():
        return None
    elif try_find.shape[0] > 1:
        raise ValueError("Multiple entries in `df` match to the key!")
    return try_find[value_col].item()


short_name2acc: dict = dict(zip(pfam_map["short_name"], pfam_map["accession"]))

data = pl.read_csv(d, separator="\t", null_values="NA").select(
    ["ProteinId", "GroupUP"] + ANNO_COLS
)


def get_pfam_acc(pfam):
    if not pfam:
        return pfam
    if REGEXES["pfam"].match(pfam):
        return REGEXES["pfam"].findall(pfam)[0]
    return short_name2acc.get(pfam, pfam)


def entry_name_from_header(header: str) -> str:
    if not header:
        return header
    if REGEXES["uniprot"].match(header):
        return REGEXES["uniprot"].findall(header)[0]
    elif REGEXES["ncbi"].match(header):
        return REGEXES["ncbi"].findall(header)[0]
    return header


ADDED = ["PFAM_IDs", "name"]


def table(lst: list):
    temp = pl.Series(lst).value_counts()
    return dict(zip(temp[""], temp["count"]))


class Grouper:

    def __init__(self, map_df: pl.DataFrame, data: pl.DataFrame) -> None:
        self.map: pl.DataFrame = map_df
        aggregated = (
            data.with_columns(
                name=pl.col("header").map_elements(
                    entry_name_from_header, return_dtype=pl.String
                )
            )
            .with_columns(pl.col(ANNO_COLS).str.split(by=";"))
            .with_columns(
                PFAM_IDs=pl.col("PFAMs").map_elements(
                    lambda x: list(map(get_pfam_acc, x)), return_dtype=pl.List(str)
                )
            )
            .group_by("GroupUP")
            .agg(pl.col("ProteinId"), pl.col(ANNO_COLS + ADDED).explode())
        )
        self.data = aggregated

    def group_count_row(self, to_count: list, source_db: str) -> pl.DataFrame:
        in_groups = {}
        for element in to_count:
            if element:
                mapped_group = get_from_df(
                    element, "Accession", "Group", self.map, ("Source", source_db)
                )
                if mapped_group:
                    in_groups[mapped_group] = in_groups.get(mapped_group, 0) + 1
        if in_groups:
            return (
                pl.DataFrame(in_groups)
                .melt(variable_name="Group")
                .sort(pl.col("value"))[-1]
            )
        return pl.DataFrame()

    def conclude_group(self, row: dict) -> str:
        """From a polars row in dictionary form, conclude what Group the
        row's entry should be assigned to based on the Group mapping provided in `map_df`
        """

        pairs = [  # Tuple of (column in row to pull elements from, column in map_df to filter on)
            ("PFAM_IDs", "PFAM"),
            ("interpro_accession", "INTERPRO"),
            ("PANTHER", "PANTHER"),
        ]
        top: pl.DataFrame = pl.DataFrame(
            {"Group": [], "value": []}, schema={"Group": pl.String, "value": pl.Int64}
        )
        for p in pairs:
            if row[p[0]] != [None]:
                found_group = self.group_count_row(row[p[0]], p[1])
                if not found_group.is_empty():
                    top = top.vstack(found_group)
        if top.is_empty():
            return "NA"
        winner = top.group_by("Group").sum().sort(pl.col("value"))[-1]
        return winner["Group"].item()

    def assign_groups(self) -> pl.DataFrame:
        groups = []
        for row in self.data.rows(named=True):
            g = self.conclude_group(row)
            groups.append(g)
        return self.data.select("ProteinId").with_columns(Group=pl.Series(groups))


# TODO: Add in the queries for JFTs and checks for eggNOG

g = Grouper(data=data, map_df=toxin_map)
grouped = g.assign_groups()

# grouped["Group"]
# grouped.filter(pl.col("Group") != "NA")
