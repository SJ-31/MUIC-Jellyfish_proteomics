#!/bin/bash

grep ">" "$1" | \
    sed 's/^>//' | \
    awk '{name=$1; $1=""}
    {
    if (name ~ /\|/)
        { split(name, ids, "|"); printf "%s\t%s\n",ids[2],$0 }
    else
        { printf "%s\t%s\n",name,$0}
     }'
