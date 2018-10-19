#!/bin/bash

stat="$1"
shift
prod=$1
shift
lnr=$1
shift
fspec=$1
shift
filter="$*"

 watch -tn 15 finfo-filter-dual.pl \"$stat\" $prod $lnr $fspec "$filter"
