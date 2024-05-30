import glob
import tempfile
import re
import pandas as pd
from subprocess import run
import os


def cleanGO(go_string):
    if not isinstance(go_string, str):
        raise ValueError(f"go_string must be a string! Received: {type(go_string)}")
    return ",".join(list(filter(lambda x: x != "NA", re.split(";", go_string))))


class Ontologizer:
    def __init__(self, go_df: pd.DataFrame, executable: str, go_path: str):
        go_df = go_df[~go_df["GO_IDs"].isna()]
        go_df["GO"] = go_df["GO_IDs"].apply(cleanGO)
        del go_df["GO_IDs"]
        go_df = go_df[go_df["GO"] != ""]
        mlist: list = ["GoStat IDs Format Version 1.0\n"]
        for id, gos in zip(go_df["ProteinId"], go_df["GO"]):
            mlist.append(f"{id}\t{gos}\n")
        self.go_path = go_path
        self.mapping: str = "".join(mlist)
        self.universe = go_df["ProteinId"]
        self.executable = executable

    def enrich(self, group_name: str, group_members: list) -> pd.DataFrame:
        cmd = [
            "java",
            "-jar",
            self.executable,
            "-g",
            self.go_path,
            "-a",
            "protein_mapping.ids",
            "-p",
            "universe.txt",
            "-s",
            f"{group_name}.txt",
        ]
        for flag, val in self.params.items():
            cmd.append(flag)
            if val != None:
                cmd.append(val)
        with open(f"{group_name}.txt", "w") as w:
            w.write("\n".join(group_members))
        run(" ".join(cmd), shell=True)
        try:
            read = pd.read_csv(glob.glob(f"table-{group_name}-*.txt")[0], sep="\t")
            return read
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    def runAll(self, groups: dict[str, list], params: dict):
        self.params = params
        results = {}
        with tempfile.TemporaryDirectory() as tmpdirname:
            pop = os.getcwd()
            os.chdir(tmpdirname)
            try:
                with open("protein_mapping.ids", "w") as w:
                    w.write(self.mapping)
                self.universe.to_csv("universe.txt", header=False, index=False)
                for group, group_members in groups.items():
                    results[group] = self.enrich(group, group_members)
            finally:
                os.chdir(pop)
        return results
