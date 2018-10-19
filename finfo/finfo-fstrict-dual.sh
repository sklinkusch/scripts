#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h $SOURCE ]; do
 DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
 SOURCE="$(readlink "$SOURCE")"
 [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

stat="$1"
prod=$2
lnr=$3

if [[ $# -eq 3 ]]; then
 watch -tn 15 $DIR/finfo-fstrict-dual.pl \"$stat\" $prod $lnr
elif [[ $# -eq 3 ]]; then
 watch -tn 15 $DIR/finfo-fstrict-dual.pl \"$stat\" $prod
fi
