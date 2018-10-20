#!/bin/bash

if [ $# -le 2 ]; then
 echo "Usage: twitter.sh <browser> <modus> <urls> ..."
 echo "<browser> = elinks or <browser> = w3m"
 echo "<modus> = s, S, m, M"
 echo "s: twitter.com + news.google.de"
 echo "m: mobile.twitter.com + news.google.de"
 echo "S: twitter.com"
 echo "M: mobile.twitter.com"
 echo "<urls> = twitter accounts (without @), separated by spaces"
 exit 0
fi
browser=$1
shift
modus=$1
shift

if [ $modus = 's' -o $modus = 'm' ]; then
 x=1
else
 x=0
fi
if [ $modus = 's' -o $modus = 'm' ]; then
 declare urls[0]=`echo "http://news.google.com/news/?ned=de&hl=de"`
fi
while [ $# -gt 0 ]; do
  b=`echo "$1"`
  if [ $modus = 'm' -o $modus = 'M' ]; then
    declare urls[$x]=`echo "http://mobile.twitter.com/${b}"`
  elif [ $modus = 's' -o $modus = 'S' ]; then
    declare urls[$x]=`echo "http://twitter.com/${b}"`
  else
    echo "Unknown modus: $modus\n"
    exit 2
  fi
  shift
  x=`expr $x + 1`
done
  outf=`perl -le 'print join " ",@ARGV' "${urls[@]}"`
  if [ $browser = 'w3m' ]; then
    options='-N'
  elif [ $browser = 'elinks' ]; then
    options=''
  else
    echo "Unknown browser: $browser\n"
    exit 3
  fi
$browser $options $outf
