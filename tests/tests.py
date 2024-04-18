#!/usr/bin/env python
import requests

entry = "https://rest.kegg.jp/list"
s = "hsa_M00001"
response = requests.get(f"{entry}/module/hsa") # Getting organism-specific
# modules seems impossible
print(response)
print(response.text)
