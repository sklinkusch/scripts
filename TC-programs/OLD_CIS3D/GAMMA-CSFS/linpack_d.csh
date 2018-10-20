#!/bin/csh
#
set echo
#
cp linpack_d.H /$HOME/include
#
g++ -c -g -I/$HOME/include linpack_d.C >& compiler.out
if ( $status != 0 ) then
  echo "Errors compiling linpack_d.C"
  exit
endif
rm compiler.out
#
mv linpack_d.o ~/libcpp/$ARCH/linpack_d.o
#
echo "A new version of linpack_d has been created."
