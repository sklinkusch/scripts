#!/bin/bash

stat="$1"
shift
mod=$1
shift
prod=$1
shift
lnr=$1
shift
fspec=$1
shift
filter="$*"

 watch -tn 15 finfo-fstrict-filter.pl \"$stat\" $mod $prod $lnr $fspec "$filter"
