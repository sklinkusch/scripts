#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h $SOURCE ]; do
 DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
 SOURCE="$(readlink "$SOURCE")"
 [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

stat="$1"
mod=$2
prod=$3
lnr=$4

if [[ $# -eq 4 ]]; then
 watch -tn 15 $DIR/finfo.pl \"$stat\" $mod $prod $lnr
elif [[ $# -eq 3 ]]; then
 watch -tn 15 $DIR/finfo.pl \"$stat\" $mod $prod
fi
