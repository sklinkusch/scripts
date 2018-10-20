#!/bin/bash

infile=`echo "$1"`
nros=`echo "$2"`

for (( i = 0; i <= $nros; i++ )); do
    if [ $i -ne 0 ]; then
	grep "STATE DIPOLE.*BOHR" "$infile" | sed -n "${i}p" | awk '{print $4, $5, $6}' 
    else
	grep "GROUND STATE (SCF) DIPOLE" "$infile" | awk '{print $5/2.5417464, $6/2.5417464, $7/2.5417464}' 
    fi
    beg=`expr $i \* $nros + 1`
    end=`expr $i \* $nros + $nros - $i`
    grep "TRANSITION DIPOLE.*BOHR" "$infile" | sed -n "$beg, ${end}p" | awk '{print $4, $5, $6}' 
done
