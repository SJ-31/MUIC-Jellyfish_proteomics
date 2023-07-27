#!/bin/bash

awk '{
if ($0 ~ />/) print $1
else print
}' "$1"
