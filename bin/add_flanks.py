#!/usr/bin/env python

import sys
modify = sys.argv[1]

with open(modify, "r") as t:
    lines = [line.split("\t") for line in t.readlines()]

for line in lines[1:]:
    current = line[103]
    line[103] = current[0] + "." + current[1:-1] + "." + current[-1]

output = ["\t".join(line) for line in lines]
with open(modify, "w") as t:
    t.write(''.join(output))
