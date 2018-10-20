#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++
F77        = gfortran
cc         = gcc
#parCC      = /home/jctremblay/bin/Libraries/mpich-3.0.4/src/env/mpicxx

## Fast flags
CFLAGS     =     -O3 -c -Wall -m64
LINK_FLAGS =     -O3 -Wall -m64 #-gfortran
LAPACK_LIBS    = -llapack -lblas
ACML_LIBS  =     -I/home/klinkusch/bin/Libraries/acml/gfortran64/include -L//home/klinkusch/bin/Libraries/acml/gfortran64/lib -L/home/klinkusch/binLibraries/acml/gfortran64_mp/lib/libacml.so -lacml -Wl,-rpath=/home/klinkusch/binLibraries/acml/gfortran64/lib/ -lgfortran -lpthread 
SCALAPACK_LIBS = -L/usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_blas95_lp64.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_lapack95_lp64.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_cdft_core.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_intel_lp64.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_sequential.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_core.a /usr/local/sw/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm 
MPI_LIBS       = -L/home/jctremblay/bin/Libraries/mpich-install/lib -lmpich -lopa -lmpl  -lrt -lpthread
OPENMP_LIBS    = -fopenmp 

## Debuging flags
#CFLAGS     =     -pg  -g   -c -Wall
#LINK_FLAGS =     -d  -pg  -g   -Wall -llapack -lblas  /usr/lib/libf2c.a

## Opteron-Cluster
#CFLAGS     =     -O3 -c -Wall
#LINK_FLAGS =     -O3 -Wall -static -L/home/software/arch64/acml/gnu64/lib/ -lacml -lg2c

OBJ_CISD  = cis-d.o ops_dcorr.o ops_cis.o ops_rhf.o ops_io.o ops_mat.o ops_crs.o
OBJ_CORG  = cis-d_orig.o ops_dcorr.o ops_cis.o ops_rhf.o ops_io.o ops_mat.o ops_crs.o
OBJ_INTS  = motwoelint.o ops_io.o ops_cis.o ops_crs.o ops_rhf.o
OBJ_RDWV  = readwav.o ops_io.o 
#OBJ_RDMT  = readmat.o ops_io.o read_gamess.o
OBJ_RDBN  = readbin.o ops_io.o
OBJ_RDSY  = readsys.o ops_io.o

all:    cis-d rgw readwav readlog wrtmat readbin readint readsys ints 
#readmat

cis-d:	$(OBJ_CISD)
	$(CC) -o ./bin/cis-d $(OBJ_CISD)  $(LINK_FLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS)

corg:	$(OBJ_CORG)
	$(CC) -o ./bin/cis-d-orig $(OBJ_CORG) $(LINK_FLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS)

ints:	$(OBJ_INTS)
	$(CC) -o ./bin/motwoelint $(OBJ_INTS) $(LINK_FLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS)

rgw:
	$(CC) -O3 -Wall -o ./bin/rgw rgw.cc

readwav:$(OBJ_RDWV)
	$(CC) -O3 -Wall -o ./bin/readwav $(OBJ_RDWV) $(LINK_FLAGS)

#readmat:$(OBJ_RDMT)
#	$(CC) -O3 -Wall -o ./bin/readmat $(OBJ_RDMT) $(LINK_FLAGS)

readlog:
	$(CC) -O3 -Wall -o ./bin/readlog readlog.cc

wrtmat:
	$(CC) -O3 -Wall -o ./bin/writemat writemat.cc

readbin:$(OBJ_RDBN)
	$(CC) -O3 -Wall -o ./bin/readbin $(OBJ_RDBN) $(LINK_FLAGS)

readint:
	$(CC) -O3 -Wall -o ./bin/readint readf08.cc

readsys:$(OBJ_RDSY)
	$(CC) -O3 -Wall -o ./bin/readsys $(OBJ_RDSY) $(LINK_FLAGS)

clean:	
	rm -f *.o *.cc~

call:	clean all

.cc.o:	
	$(CC) $(CFLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS) $<

.f.o:
	$(F77) -fdefault-integer-8 -frecord-marker=8 -c   $<
#.f.o:
#	$(F77) -c   $<

.c.o:	
	$(cc) -c -O3 -Wall $<
