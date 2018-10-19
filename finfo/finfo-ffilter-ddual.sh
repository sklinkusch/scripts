#!/bin/bash

stat="$1"
shift
mod=$1
shift
prod=$1
shift
lnr=$1
shift
faspec=$1
shift
fbspec=$1
shift
filter="$*"

 watch -tn 15 finfo-ffilter-ddual.pl \"$stat\" $mod $prod $lnr $faspec $fbspec "$filter"
