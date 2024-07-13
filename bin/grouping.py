import polars as pl
import tomllib
import re

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


def get_pfam_acc(pfam, mapping):
    if not pfam:
        return pfam
    if REGEXES["pfam"].match(pfam):
        return REGEXES["pfam"].findall(pfam)[0]
    return mapping.get(pfam, pfam)


def entry_name_from_header(header: str) -> str:
    if not header:
        return header
    if REGEXES["uniprot"].match(header):
        return REGEXES["uniprot"].findall(header)[0]
    elif REGEXES["ncbi"].match(header):
        return REGEXES["ncbi"].findall(header)[0]
    return header


ADDED = ["PFAM_IDs", "name"]


class Grouper:

    def __init__(
        self, map_df: pl.DataFrame, data: pl.DataFrame, pfam_map_path: str
    ) -> None:
        pfam_map = pl.read_csv(pfam_map_path, separator="\t")
        short_name2acc: dict = dict(zip(pfam_map["short_name"], pfam_map["accession"]))

        self.map: dict = dict(zip(map_df["Accession"], map_df["Group"]))
        aggregated = (
            data.with_columns(
                name=pl.col("header").map_elements(
                    entry_name_from_header, return_dtype=pl.String
                )
            )
            .with_columns(pl.col(ANNO_COLS).str.split(by=";"))
            .with_columns(
                PFAM_IDs=pl.col("PFAMs").map_elements(
                    lambda x: list(map(lambda y: get_pfam_acc(y, short_name2acc), x)),
                    return_dtype=pl.List(str),
                )
            )
            .group_by("GroupUP")
            .agg(pl.col("ProteinId"), pl.col(ANNO_COLS + ADDED).explode())
        )
        self.cols_to_count_groups = ["PFAM_IDs", "interpro_accession", "PANTHER"]
        # columns containing accessions that will be counted to determine groups
        self.cols_to_sum = self.cols_to_count_groups.copy()
        # names of columns containing the groups that will be considered when summarizing the group for an entry in  the `assign_groups` function
        self.data = aggregated

    def group_count_row(self, to_count: list) -> tuple:
        group_tracker = {}
        for element in to_count:
            if element:
                mapped_group = self.map.get(element)
                if mapped_group:
                    group_tracker[mapped_group] = group_tracker.get(mapped_group, 0) + 1
        if group_tracker:
            top = sorted(group_tracker.items(), key=lambda x: x[1])[-1]
            return {"group": top[0], "count": top[1]}
        return {"group": "NA", "count": 0}

    def assign_groups(self) -> pl.DataFrame:
        """From a polars row in dictionary form, conclude what Group the
        row's entry should be assigned to based on the Group mapping provided in `map_df`
        """
        concluded = self.data.with_columns(
            pl.col(self.cols_to_count_groups).map_elements(
                self.group_count_row,
                return_dtype=pl.Struct({"group": pl.String, "count": pl.Int64}),
            ),
        )
        concluded = concluded.select(
            "GroupUP",
            pl.struct(self.cols_to_sum)
            .map_elements(find_max_struct, return_dtype=pl.String)
            .alias("Group"),
        )
        return concluded.filter(pl.col("Group") != "")


def find_max_struct(x):
    return max(x.values(), key=lambda x: x["count"])["group"]


def mapping_from_definition(dct: dict[str, list]) -> dict:
    """Maps the elements of the values of dct to their keys"""
    return {item: key for key, value in dct.items() for item in value}


def find_in_headers(headers, mapping: dict) -> dict:
    tracker = {}

    def find_one(header: str):
        for k, v in mapping.items():
            if re.match(k.lower(), header.lower()):
                tracker[v] = tracker.get(v, 0) + 1

    [find_one(q) for q in headers]
    if not tracker:
        return {"group": "NA", "count": 0}

    top = sorted(tracker.items(), key=lambda x: x[1])[-1]
    return {"group": top[0], "count": top[1]}


class ToxinGrouper(Grouper):

    def __init__(
        self,
        map_df: pl.DataFrame,
        data: pl.DataFrame,
        toxin_header_map: dict,
        pfam_map_path: str,
    ) -> None:
        super().__init__(map_df, data, pfam_map_path)
        self.cols_to_sum.append("group_from_header")
        self.toxin_header_map = mapping_from_definition(toxin_header_map)

    def assign_groups(self) -> pl.DataFrame:
        self.data = self.data.with_columns(
            group_from_header=pl.col("name").map_elements(
                lambda x: find_in_headers(x, self.toxin_header_map),
                return_dtype=pl.Struct({"group": pl.String, "count": pl.Int64}),
            )
        )
        return super().assign_groups()


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-p", "--pfam_map_path")
    parser.add_argument("-t", "--toxin_map_path")
    parser.add_argument("-a", "--all_map_path")
    parser.add_argument("-m", "--mode")
    args = vars(parser.parse_args())
    return args


def main(args):
    data = pl.read_csv(args["input"], separator="\t", null_values="NA").select(
        ["ProteinId", "GroupUP"] + ANNO_COLS
    )
    toxin_map = pl.read_csv(args["toxin_map_path"], separator="\t")
    with open(args["all_map_path"], "rb") as t:
        maps = tomllib.load(t)
    if args["mode"] == "toxin":
        G = ToxinGrouper(
            toxin_map, data, maps["header_names"]["toxins"], args["pfam_map_path"]
        )
        result = G.assign_groups()
        result.write_csv(args["output"], separator="\t")
    elif args["mode"] == "cog":
        G = Grouper()


if __name__ == "__main__":
    args = parse_args()
    main(args)
