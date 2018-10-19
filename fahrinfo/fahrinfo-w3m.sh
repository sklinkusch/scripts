#!/bin/bash
path='/home/stefan/bin/fahrinfo-elinks.dat'
if [ $# -lt 6 -o $# -gt 9 ]; then
 echo "Gebrauch: fahrnfo-w3m.sh \"<Haltestelle>\" <#Fahrten> <dep/arr> <code> <datecode> <begindate (datecode=1)> <enddate (datecode=1)> <timecode> <Zeile (optional)>"
 echo "dep: Abfahrtsplan, arr: Ankunftsplan"
 echo "code: Summe von Zahlen: DB-Fernverkehr (64), Regionalverkehr (32), S-Bahn (16), U-Bahn (8), Tram (4), Bus (2), Faehre (1)"
 echo "datecode: 1: Zeitraum, 0: heute; begindate/enddate im Format TT.MM.YY, timecode im Format HH:MM oder 0 (aktuelle Zeit)"
 exit
fi
if [ $5 -eq 1 ]; then
    datum="period"
    begindate=$6
    enddate=$7
    timecode=$8
    if [[ $9 ]]; then
        line=$9
    else
        line='0'
    fi
else
    datum="today"
    timecode=$6
    if [[ $7 ]]; then
        line=$7
    else
        line='0'
    fi
fi
haltst="$1"
nroh=`grep -i "$haltst" $path | wc -l`
if [ $nroh -eq 0 ]; then
 echo "Keine Haltestelle gefunden"
 exit
elif [ $nroh -eq 1 ]; then
 haltestelle=`grep -i "$haltst" $path | cut -d '	' -f 2`
elif [ $nroh -gt 1 -a $nroh -le 50 ]; then
 echo "$nroh Haltestellen gefunden: "
 grep -i "$haltst" $path | cut -d '	' -f 1 | nl
 if [ $line != '0' ]; then
     haltestelle=`grep -i "$haltst" $path | sed -n "${line}p" | cut -d '	' -f 2`
 else
  exit
 fi
elif [ $nroh -gt 50 ]; then
 echo "$nroh Haltestellen gefunden! Präzisieren Sie Ihre Eingabe!"
 exit
else
 echo "Fehler"
 exit
fi
add="00"
haltestelle="${haltestelle:0:1}00${haltestelle:1}"
nroj=$2
typus=$3
filtre=$4
zero=0
one=1
if [ $filtre -le 0 -o $filtre -gt 127 ]; then
 filtre=127
fi
if [ $filtre -ge 64 ]; then
 db=1
 filtre=`expr $filtre - 64`
else
 db=0
fi
if [ $filtre -ge 32 ]; then
 rb=1
 filtre=`expr $filtre - 32`
else
 rb=0
fi
if [ $filtre -ge 16 ]; then
 sb=1
 filtre=`expr $filtre - 16`
else
 sb=0
fi
if [ $filtre -ge 8 ]; then
 ub=1
 filtre=`expr $filtre - 8`
else
 ub=0
fi
if [ $filtre -ge 4 ]; then
 tram=1
 filtre=`expr $filtre - 4`
else
 tram=0
fi
if [ $filtre -ge 2 ]; then
 bus=1
 filtre=`expr $filtre - 2`
else
 bus=0
fi
if [ $filtre -eq 1 ]; then
 ferry=1
else
 ferry=0
fi

filter=$sb$ub$tram$bus$ferry$db$rb$zero$one
if [ $timecode = '0' -a $datum = "today" ]; then
    url=`echo "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestelle&boardType=$typus&maxJourneys=$nroj&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&"`
elif [ $datum = "today" ]; then
    url=`echo "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestelle&boardType=$typus&time=$timecode&maxJourneys=$nroj&selectDate=today&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&"`
elif [ $timecode = '0' ]; then
    url=`echo "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestelle&boardType=$typus&maxJourneys=$nroj&selectDate=period&dateBegin=$begindate&dateEnd=$enddate&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&"`
else
    url=`echo "http://fahrinfo.vbb.de/bin/stboard.exe/dn?input=$haltestelle&boardType=$typus&time=$timecode&maxJourneys=$nroj&selectDate=period&dateBegin=$begindate&dateEnd=$enddate&productsFilter=$filter&start=yes&pageViewMode=PRINT&dirINPUT=&"`
fi

 w3m "$url"

