#!/usr/bin/env python

import pandas as pd
import scipy.cluster.hierarchy as scipy


def findNewGroups(df: pd.DataFrame) -> None:
    groupings = scipy.DisjointSet()
    gs = set()
    id_to_group = {}
    # Establish groups and find roots
    for pid, t in zip(df["ProteinId"], df["ProteinGroupId"]):
        splits = t.split(";")
        first = splits[0]
        gs.add(first)
        if first not in groupings:
            groupings.add(first)
        id_to_group[pid] = first
        if len(splits) > 1:
            for s in splits[1:]:
                if s not in groupings:
                    groupings.add(s)
                groupings.merge(first, s)

    count = 0
    # Map proteins to new roots
    group_key: dict = {}
    new_groups: pd.Series = []
    for k, v in id_to_group.items():
        root = groupings.__getitem__(v)
        if root in group_key:
            group = group_key.get(root)
        else:
            group = f"G{count}"
            group_key[root] = group
            count += 1
        new_groups.append(group)
    new_groups = pd.Series(new_groups)
    df["Group"] = new_groups
    return df
