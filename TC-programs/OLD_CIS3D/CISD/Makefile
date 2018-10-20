#********************************************************************************
#* CISD Program suite                                                           *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                      Tillmann Klamroth,Pascal Krause  2005   *
#********************************************************************************/ 

CC         = g++
F77        = gfortran-4.1
cc         = gcc

## Opteron-Cluster
CFLAGS     =     -O3 -c -Wall
LINK_FLAGS =     -O3 -Wall -llapack -lblas -lgfortran


OBJ_CISD  = cisd.o ops_io.o ops_mat.o ops_cisd.o ops_rhf.o ops_slater_condon.o ops_mk_csf.o
OBJ_TDCISD = td_cisd.o ops_io.o ops_mat.o ops_las.o ops_tdcisd.o

all:    cisd tdcisd

cisd:	$(OBJ_CISD)
	$(CC) -o ../bin/cisd   $(OBJ_CISD) $(LINK_FLAGS) 

tdcisd:	$(OBJ_TDCISD)
	$(CC) -o ../bin/tdcisd $(OBJ_TDCISD) $(LINK_FLAGS)

clean:	
	rm -f *.o *.cc~

call:	clean all


.cc.o:	
	$(CC) $(CFLAGS) $<

