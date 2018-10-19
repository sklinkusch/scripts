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
amodus=$1
shift
aprod=$1
shift
alnr=$1
shift
bstat="$1"
shift
bmodus=$1
shift
bprod=$1
shift
blnr=$1

watch -tn 15 $DIR/finfo-fstrict-duo.pl \"$astat\" $amodus $aprod $alnr \"$bstat\" $bmodus $bprod $blnr 
