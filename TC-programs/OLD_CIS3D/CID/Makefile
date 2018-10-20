#********************************************************************************
#* CISD Program suite                                                           *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                      Tillmann Klamroth,Pascal Krause  2005   *
#********************************************************************************/ 

CC         = g++
F77        = g77
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


OBJ_CID  = cid.o ops_io.o ops_mat.o ops_cid.o ops_rhf.o ops_slater_condon.o ops_mk_csf.o

all:    cid 

cid:	$(OBJ_CID)
	$(CC) -o cid   $(OBJ_CID) $(LINK_FLAGS) 

clean:	
	rm -f *.o *.cc~

call:	clean all


.cc.o:	
	$(CC) $(CFLAGS) $<

