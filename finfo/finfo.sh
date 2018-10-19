#!/bin/bash

stat="$1"
mod=$2
prod=$3
lnr=$4

if [[ $# -eq 4 ]]; then
 watch -tn 15 finfo.pl \"$stat\" $mod $prod $lnr
elif [[ $# -eq 3 ]]; then
 watch -tn 15 finfo.pl \"$stat\" $mod $prod
fi
