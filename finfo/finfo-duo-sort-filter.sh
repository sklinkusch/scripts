#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h $SOURCE ]; do
 DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
 SOURCE="$(readlink "$SOURCE")"
 [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

astat="$1"
shift
bstat="$1"
shift
mod=$1
shift
prod=$1
shift
alnr=$1
shift
blnr=$1
shift
fspec=$1
shift
filter="$*"

 watch -tn 15 $DIR/finfo-duo-sort-filter.pl \"$astat\" \"$bstat\" $mod $prod $alnr $blnr $fspec "$filter"
