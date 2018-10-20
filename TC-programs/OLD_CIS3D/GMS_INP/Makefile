#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++
F77        = gfortran-4.3
cc         = gcc

## Fast flags
CFLAGS     =     -O3 -c -Wall
LINK_FLAGS =     -O3 -Wall -llapack -lblas -lgfortran

## Debuging flags
#CFLAGS     =     -pg  -g   -c -Wall
#LINK_FLAGS =     -d  -pg  -g   -Wall -llapack -lblas  /usr/lib/libf2c.a

## Opteron-Cluster
#CFLAGS     =     -O3 -c -Wall
#LINK_FLAGS =     -O3 -Wall -static -L/home/software/arch64/acml/gnu64/lib/ -lacml -lg2c


OBJ_GMS   = main.o dft.o basis.o contrl.o ordz.o

all:    gms 

gms:    $(OBJ_GMS)
	$(CC) -o gms_inp  $(OBJ_GMS)  $(LINK_FLAGS)

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $<

.f.o:
	$(F77)  -c $<

.c.o:	
	$(cc) -c -O3 -Wall $<
