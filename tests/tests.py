#!/usr/bin/env ipython
import pytest
from pathlib import Path
import xml.etree.ElementTree as ET
import pandas as pd
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
test_path = "/home/sc31/Bio_SDD/MUIC_senior_project/workflow/unittests"
sys.path.append("/home/sc31/Bio_SDD/MUIC_senior_project/workflow/bin")
from make_manifest import make_manifest
from merge_quantifications import merge_dlfq,read_directlfq
from emPAI import calculate_emPAI

def test_make_manifest():
    template = pd.DataFrame({"Raw": ["one.raw", "two.raw"],
                             "indexed_mzML": ["one.mzML", "two.mzML"],
                             "mgf": ["one.mgf", "two.mgf"]})
    expected_mzml = [f"{test_path}/empty_files/three.mzML",
                     f"{test_path}/empty_files/four.mzML"]
    expected_mgf = [f"{test_path}/empty_files/six.mgf",
                    f"{test_path}/empty_files/five.mgf"]
    test = make_manifest(template, f"{test_path}/empty_files")
    assert list(test["indexed_mzML"]) == expected_mzml
    assert list(test["mgf"]) == expected_mgf

# test_make_manifest()

class MergeQuant:

    def __init__(self, file_path) -> None:
        dlfq = f'{file_path}/1-First_pass/Quantify/directlfq_prot.tsv'
        seq_mapping = f'{file_path}/Databases/all_normal_mapping.tsv'
        intersected = f'{file_path}/1-First_pass/Combined/unified_groups.tsv'
        self.dlfq = read_directlfq(dlfq)
        self.mapping = pd.read_csv(seq_mapping, sep="\t")
        self.proteins = pd.read_csv(intersected, sep="\t")
        self.merged = merge_dlfq(dlfq, self.proteins, float(1))
        self.qdf = calculate_emPAI(self.merged, self.mapping, (350, 1600))


path = "../results/ND_jellyfish"
# mq = MergeQuant(path)

dlfq = f'{path}/1-First_pass/Quantify/directlfq_prot.tsv'
# df = read_directlfq(dlfq)
#

from maxquant_wrapper import MQConfig
mq_config = MQConfig("./mq_directory/default_maxquant.xml")
db = str(Path("./mq_directory/test_db.fasta").absolute())
raws = list(Path("./mq_directory").glob("*.raw"))
mq_config.add_fasta(db)
mq_config.add_paths(raws)
mq_config.write_file("./mq_directory/all_raws.xml")
