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


OBJ_UHF   = uhf.o ops_io.o ops_mat.o ops_mat_neu.o ops_uhf.o
OBJ_CIS   = cis.o ops_io.o ops_mat.o ops_mat_neu.o ops_uhf.o ops_cis.o ops_so_corr.o ops_ro_corr.o
OBJ_TDRHF = td_rhf.o ops_io.o ops_mat.o ops_mat_neu.o ops_uhf.o ops_tdrhf.o ops_las.o 
OBJ_TDCIS = td_cis.o ops_io.o ops_mat.o ops_mat_neu.o ops_las.o ops_tdcis.o
OBJ_TDPOP = td_pop.o ops_io.o ops_uhf.o ops_mat.o ops_mat_neu.o ops_popa.o 
OBJ_MKIN  = mk_in.o  read_gamess.o
OBJ_MKIN07 = mk_in_neu.o read_gms07.o 
OBJ_MKIW  = mk_init_wav.o ops_io.o ops_mat.o ops_mat_neu.o ops_uhf.o ops_ro_corr.o ops_cis.o
OBJ_LCORB = loc_orbs.o ops_io.o ops_mat.o ops_mat_neu.o ops_uhf.o
OBJ_OCT   = cis_oct.o ops_oct.o ops_io.o ops_mat.o ops_mat_neu.o ops_tdcis.o ops_mt.o

all:    uhf cis tdrhf tdcis tdpop mkin mkin07 rgw mkiw lcorb mkef oct 

uhf:	$(OBJ_UHF)
	$(CC) -o ../uhf-bin/uhf   $(OBJ_UHF)   $(LINK_FLAGS) 

cis:	$(OBJ_CIS)
	$(CC) -o ../uhf-bin/cis   $(OBJ_CIS)   $(LINK_FLAGS) 

tdrhf:	$(OBJ_TDRHF)
	$(CC) -o ../uhf-bin/tdrhf $(OBJ_TDRHF) $(LINK_FLAGS)

tdcis:	$(OBJ_TDCIS)
	$(CC) -o ../uhf-bin/tdcis $(OBJ_TDCIS) $(LINK_FLAGS)	

tdpop:	$(OBJ_TDPOP)
	$(CC) -o ../uhf-bin/tdpop $(OBJ_TDPOP) $(LINK_FLAGS)

mkin:	$(OBJ_MKIN)
	$(CC) -o ../uhf-bin/mk_in $(OBJ_MKIN)  $(LINK_FLAGS)

mkin07:	$(OBJ_MKIN07)
	$(CC) -o ../uhf-bin/mk_in07 $(OBJ_MKIN07)  $(LINK_FLAGS)

rgw:
	$(CC) -O3 -Wall -o ../bin/rgw rgw.cc

mkiw:	$(OBJ_MKIW)
	$(CC) -o ../bin/mk_init_wav $(OBJ_MKIW)  $(LINK_FLAGS)

lcorb:	$(OBJ_LCORB)
	$(CC) -o ../bin/lcorb	$(OBJ_LCORB) $(LINK_FLAGS)

oct:	$(OBJ_OCT)
	$(CC) -o ../bin/cis_oct	$(OBJ_OCT) $(LINK_FLAGS) -lpthread

mkef:	mk_field.cc
	$(CC) -o ../bin/mk_field -O3 mk_field.cc -lfftw

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $<

.f.o:
	$(F77)  -c $<

.c.o:	
	$(cc) -c -O3 -Wall $<
