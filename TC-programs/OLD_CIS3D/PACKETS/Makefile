#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++
cc         = gcc

## Fast flags
CFLAGS     =     -O3 -c -Wall
LINK_FLAGS =     -O3 -Wall -llapack -lblas -lgfortran



OBJ_MKIW  = mk_init_wav.o ops_io.o ops_mat.o ops_rhf.o ops_ro_corr.o ops_cis.o
OBJ_OCT   = cis_oct.o ops_oct.o ops_io.o ops_mat.o ops_tdcis.o ops_mt.o

all:     mkiw  oct mkef 


mkiw:	$(OBJ_MKIW)
	$(CC) -o mk_init_wav $(OBJ_MKIW)  $(LINK_FLAGS)


oct:	$(OBJ_OCT)
	$(CC) -o cis_oct	$(OBJ_OCT) $(LINK_FLAGS) -lpthread

mkef:	mk_field.cc
	$(CC) -o mk_field -O3 mk_field.cc -lfftw

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $<

.c.o:	
	$(cc) -c -O3 -Wall $<
