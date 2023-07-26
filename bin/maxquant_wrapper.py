#!/usr/bin/python3
"""Wrapper for specifying input files to maxqant"""
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
from pathlib import Path

parser = ArgumentParser(prog="MaxQuant Wrapper",
                        description="Specify input files to MaxQuant config")
parser.add_argument("raw_file")
parser.add_argument("config_template")
args = parser.parse_args()
raw = Path(args.raw_file)

config = ET.parse(args.config_template)
root = config.getroot()
root.set("xmlns:xsd", "http://www.w3.org/2001/XMLSchema")
root.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
print(root)
filepaths = root.find("filePaths")

path = filepaths.find("string")
if path is not None:
    filepaths.remove(path)
#   Now add all paths
sub = ET.SubElement(filepaths, "string")
sub.text = str(raw.absolute())
config.write("mqconfig.xml",
             encoding="utf-8",
             method="xml",
             xml_declaration=True)
