#!/bin/bash

stat="$1"
shift
mod=$1
shift
prod=$1
shift
lnr=$1

 watch -tn 15 finfo-strict.pl \"$stat\" $mod $prod $lnr
