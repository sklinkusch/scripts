#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h $SOURCE ]; do
 DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
 SOURCE="$(readlink "$SOURCE")"
 [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

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

 watch -tn 15 $DIR/finfo-fstrict-filter.pl \"$stat\" $mod $prod $lnr $fspec "$filter"
