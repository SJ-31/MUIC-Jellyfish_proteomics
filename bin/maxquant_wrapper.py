#!/usr/bin/python3
"""Wrapper for specifying input files to maxqant"""
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from pathlib import Path


class MQConfig:

    def __init__(self, root_path) -> None:
        self.config = ET.parse(root_path)
        self.root = self.config.getroot()
        self.root.set("xmlns:xsd", "http://www.w3.org/2001/XMLSchema")
        self.root.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")

    def add_fasta(self, fasta):
        for f in self.root.find("fastaFiles").iter():
            if f.tag == "fastaFilePath":
                f.text = fasta

    def add_paths(self, files,
                  exp_text, frac_text,
                  index_text, ref_text,
                  ptm_text):
        filepaths_xml = self.root.find("filePaths")
        experiment_xml = self.root.find("experiments")
        fractions_xml = self.root.find("fractions")
        indices_xml = self.root.find("paramGroupIndices")
        reference_channel_xml = self.root.find("referenceChannel")
        ptms_xml = self.root.find("ptms")
        old_path = filepaths_xml.find("string")
        if old_path is not None:
            filepaths_xml.remove(old_path)
        #   Now add all paths
        for file in files:
            sub = ET.SubElement(filepaths_xml, "string")
            sub.text = str(file.absolute())
            frac = ET.SubElement(fractions_xml, "short")
            frac.text = "32767"
            exp = ET.SubElement(experiment_xml, "string")
            exp.text = exp_text
            index = ET.SubElement(indices_xml, "index")
            index.text = index_text
            ref = ET.SubElement(reference_channel_xml, "string")
            ref.text = ref_text
            ptm = ET.SubElement(ptms_xml, "boolean")
            ptm.text = ptm_text

    def write_file(self, file_name):
        self.config.write(file_name,
                          encoding="utf-8",
                          method="xml",
                          xml_declaration=True)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser = ArgumentParser(
        prog="MaxQuant Wrapper",
        description="Specify input files to MaxQuant config")
    parser.add_argument("-r,", "--raw_file_path")
    parser.add_argument("-c", "--config_template")
    parser.add_argument("-o", "--output")
    parser.add_argument("-d", "--database_path")
    parser.add_argument("--ptms")  # Default: False
    parser.add_argument("--fraction")  # Default: 32767
    parser.add_argument("--index")  # Default: 0
    parser.add_argument("--experiment")  # Default: None ("")
    parser.add_argument("--reference_channel")  # Default: None ("")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == '__main__':
    args = parse_args()
    raws = Path(args["raw_file_path"]).glob("*.raw")
    db = str(Path(args["database_path"]).absolute())
    config = MQConfig(args["config_template"])
    config.add_fasta(db)
    config.add_paths(files=raws, frac_text=args["fraction"],
                     index_text=args["index"], exp_text=args["experiment"],
                     ref_text=args["reference_channel"],
                     ptm_text=args["ptms"])
    config.write_file(args["output"])
