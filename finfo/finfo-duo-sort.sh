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
modus=$1
shift
prod=$1
shift
alnr=$1
shift
blnr=$1

watch -tn 15 $DIR/finfo-duo-sort.pl \"$astat\" \"$bstat\" $modus $prod $alnr $blnr 
