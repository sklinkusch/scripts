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


OBJ_RHF   = rhf.o ops_io.o ops_mat.o ops_rhf.o
OBJ_CIS   = cis.o ops_io.o ops_mat.o ops_rhf.o ops_cis.o ops_so_corr.o ops_ro_corr.o
OBJ_TDRHF = td_rhf.o ops_io.o ops_mat.o ops_rhf.o ops_tdrhf.o ops_las.o 
OBJ_TDCIS = td_cis.o ops_io.o ops_mat.o ops_las.o ops_tdcis.o
OBJ_TDPOP = td_pop.o ops_io.o ops_rhf.o ops_mat.o ops_popa.o 
OBJ_MKIN  = mk_in.o  read_gamess.o
OBJ_MKIW  = mk_init_wav.o ops_io.o ops_mat.o ops_rhf.o ops_ro_corr.o ops_cis.o
OBJ_LCORB = loc_orbs.o ops_io.o ops_mat.o ops_rhf.o
OBJ_OCT   = cis_oct.o ops_oct.o ops_io.o ops_mat.o ops_tdcis.o ops_mt.o

all:    rhf cis tdrhf tdcis tdpop mkin rgw mkiw lcorb mkef oct

rhf:	$(OBJ_RHF)
	$(CC) -o ./bin/rhf   $(OBJ_RHF)   $(LINK_FLAGS) 

cis:	$(OBJ_CIS)
	$(CC) -o ./bin/cis   $(OBJ_CIS)   $(LINK_FLAGS) 

tdrhf:	$(OBJ_TDRHF)
	$(CC) -o ./bin/tdrhf $(OBJ_TDRHF) $(LINK_FLAGS)

tdcis:	$(OBJ_TDCIS)
	$(CC) -o ./bin/tdcis $(OBJ_TDCIS) $(LINK_FLAGS)	

tdpop:	$(OBJ_TDPOP)
	$(CC) -o ./bin/tdpop $(OBJ_TDPOP) $(LINK_FLAGS)

mkin:	$(OBJ_MKIN)
	$(CC) -o ./bin/mk_in $(OBJ_MKIN)  $(LINK_FLAGS)

rgw:
	$(CC) -O3 -Wall -o ./bin/rgw rgw.cc

mkiw:	$(OBJ_MKIW)
	$(CC) -o ./bin/mk_init_wav $(OBJ_MKIW)  $(LINK_FLAGS)

lcorb:	$(OBJ_LCORB)
	$(CC) -o ./bin/lcorb	$(OBJ_LCORB) $(LINK_FLAGS)

oct:	$(OBJ_OCT)
	$(CC) -o ./bin/cis_oct	$(OBJ_OCT) $(LINK_FLAGS) -lpthread

mkef:	mk_field.cc
	$(CC) -o ./bin/mk_field -O3 mk_field.cc -lfftw

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $<

.f.o:
	$(F77)  -c $<

.c.o:	
	$(cc) -c -O3 -Wall $<
