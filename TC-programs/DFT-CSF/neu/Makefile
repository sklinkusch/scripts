#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++
F77        = gfortran-4.5
cc         = gcc

## Fast flags
CFLAGS     =     -O3 -c -Wall -pg -g
LINK_FLAGS =     -O3 -Wall -llapack -lblas -lgfortran -pg -g
OPENMP_FLAG=     -fopenmp

## Debuging flags
#CFLAGS     =     -pg  -g   -c -Wall
#LINK_FLAGS =     -d  -pg  -g   -Wall -llapack -lblas  /usr/lib/libf2c.a

## Opteron-Cluster
#CFLAGS     =     -O3 -c -Wall
#LINK_FLAGS =     -O3 -Wall -static -L/home/software/arch64/acml/gnu64/lib/ -lacml -lg2c


OBJ_CIS   = istochdft.o ops_io.o ops_mat.o ops_rhf.o ops_cis.o ops_crs.o

OBJ_RHO   = irhodft.o ops_io.o ops_mat.o ops_rhf.o ops_cis.o ops_crs.o

OBJ_IRX   = iirxdft.o ops_io.o

all:    cis rho irx

cis:	$(OBJ_CIS)
	$(CC) -o ../bin/istochdft2   $(OBJ_CIS)   $(LINK_FLAGS) $(OPENMP_FLAG)

rho:	$(OBJ_RHO)
	$(CC) -o ../bin/irhodft2 $(OBJ_RHO) $(LINK_FLAGS) $(OPENMP_FLAG)	

irx:	$(OBJ_IRX)
	$(CC) -o ../bin/iirxdft2 $(OBJ_IRX) $(LINK_FLAGS) $(OPENMP_FLAG)

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $(OPENMP_FLAG) $<

.f.o:
	$(F77)  -c $<

.c.o:	
	$(cc) -c -O3 -Wall $<
