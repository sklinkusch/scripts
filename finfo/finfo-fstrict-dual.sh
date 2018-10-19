#!/bin/bash

stat="$1"
prod=$2
lnr=$3

if [[ $# -eq 3 ]]; then
 watch -tn 15 finfo-fstrict-dual.pl \"$stat\" $prod $lnr
elif [[ $# -eq 3 ]]; then
 watch -tn 15 finfo-fstrict-dual.pl \"$stat\" $prod
fi
