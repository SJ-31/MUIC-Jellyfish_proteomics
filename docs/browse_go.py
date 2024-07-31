#!/usr/bin/env python
import sys
import re


def see(data, **kwargs):
    import dtale
    import pandas as pd
    import polars
    import subprocess

    if isinstance(data, pd.DataFrame) or (
        is_polars := isinstance(data, polars.dataframe.frame.DataFrame)
    ):
        if is_polars:
            data = data.to_pandas()
        d = dtale.show(data)
        subprocess.run(f"firefox --new-window {d._url}", shell=True)
    elif (
        isinstance(data, list)
        or isinstance(data, dict)
        or isinstance(data, str)
        or isinstance(data, set)
        or isinstance(data, polars.series.series.Series)
    ):
        import pprint as pp
        import pathlib
        import datetime

        if isinstance(data, polars.series.series.Series):
            data = list(data)

        tmp_see = pathlib.Path("/tmp/see_tmp")
        if not tmp_see.exists():
            tmp_see.mkdir()
        time = datetime.datetime.now().strftime("%Y-%m-%d-%M-%S")
        if kwargs.get("yaml", False):
            import yaml

            tmp_file = f"{tmp_see}/{time}_see_temp.yaml"
            with open(tmp_file, "w") as f:
                yaml.dump(data, f)
        elif kwargs.get("md", False):
            tmp_file = f"{tmp_see}/{time}_see_temp.md"
            with open(tmp_file, "w") as f:
                f.write(data)
        subprocess.Popen(["nvim-qt", tmp_file])


def md_format(browse_result: list):
    def format_entry(d: dict):
        nc = d["M"].get("n_children")
        ontology = d["M"].get("ontology")
        level = d["M"].get("level")
        id = d["I"].replace('"', "")
        return f"# {id}\n**level:** {level} | **ontology:** {ontology} | **n_children:** {nc}\n{d['M'].get("definition")}\n"

    strs = [format_entry(s) for s in browse_result]
    return "\n".join(strs)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--metadata_path", default="../data/reference/with_levels.tsv"
    )
    parser.add_argument("-g", "--go_obo", default="../data/reference/go_networkx.gml")
    parser.add_argument("-f", "--format", default="yaml")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    import readline

    sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
    import go_subset as gs

    readline.parse_and_bind("set editing-mode emacs")

    def show_help() -> str:
        print("Start your command with 'b' to perform a keyword search on the GO")
        print("Optionally specify ontology with 'o=BP|CC|MF'")
        print(
            "Whitespace on sides is automatically stripped, use '_' if you want to prepend/append spaces"
        )
        print(
            "Start your command with 's' to list the successors (children) of the specified GO term"
        )
        print(
            "Start your command with 'p' to list the parents of the specified GO term"
        )

    GO = gs.CompleteGO(metadata_path=args["metadata_path"], go_path=args["go_obo"])
    browsing = True
    while browsing:
        response = input("GO CLI browser (? for help, Q to quit) \n  ")
        if response == "?":
            show_help()
        elif response.startswith("b"):
            response = response[1:]
            ontology = re.findall("o=([A-Za-z]+)", response)
            if ontology:
                response = (
                    re.sub("o=([A-Za-z]+)", "", response).strip().replace("_", " ")
                )
                o: str = ontology[0].upper()
                print(f'Searching for "{response}"...')
                result = GO.browse(response, sort=True, ontology=o)
            else:
                result = GO.browse(response, sort=True)
            if args["format"] == "yaml":
                see(result, yaml=True)
            elif args["format"] == "md":
                see(md_format(result), md=True)
        elif response.startswith("s"):
            response = response[1:].strip()
            see(GO.nested_successors(response), yaml=True)
        elif response.startswith("p"):
            response = response[1:].strip()
            see(GO.get_parents(response), yaml=True)
        elif response.startswith("Q"):
            print("Browsing session ended")
            browsing = False
