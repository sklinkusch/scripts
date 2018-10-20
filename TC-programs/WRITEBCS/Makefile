#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++
F77        = gfortran-4.7 
cc         = gcc

## Fast flags
CFLAGS     =     -O3 -c -Wall
LINK_FLAGS =     -O3 -Wall -llapack -lblas -lgfortran
OPENMP_FLAGS =   -fopenmp

## Debuging flags
#CFLAGS     =     -pg  -g   -c -Wall
#LINK_FLAGS =     -d  -pg  -g   -Wall -llapack -lblas  /usr/lib/libf2c.a

## Opteron-Cluster
#CFLAGS     =     -O3 -c -Wall
#LINK_FLAGS =     -O3 -Wall -static -L/home/software/arch64/acml/gnu64/lib/ -lacml -lg2c

OBJ_WRTBCS= write_bcs.o ops_io.o ops_cis.o ops_rhf.o ops_mat.o read_gamess.o

all:    wrtbcs 

wrtbcs:	$(OBJ_WRTBCS)
	$(CC) -o ./bin/writebcs $(OBJ_WRTBCS) $(LINK_FLAGS) $(OPENMP_FLAGS)

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) $<

.f.o:
	$(F77) -fdefault-integer-8 -frecord-marker=8 -c   $<
#.f.o:
#	$(F77) -c   $<

.c.o:	
	$(cc) -c -O3 -Wall $<
