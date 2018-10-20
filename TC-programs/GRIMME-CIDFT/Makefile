#********************************************************************************
#* CIS3(D) Program suite                                                        *
#*                                                                              *
#* file: Makefile                                                               *
#*                                                                              *
#*                                                    Tillmann Klamroth  2004   *
#********************************************************************************/ 

CC         = g++-4.7
F77        = gfortran-4.7 
cc         = gcc-4.7
#parCC      = /home/jctremblay/bin/Libraries/mpich-3.0.4/src/env/mpicxx

## Fast flags
CFLAGS     =     -O3 -c -Wall -m64
LINK_FLAGS =     -O3 -Wall -lgfortran -m64
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

OBJ_GRIM  = grimme.o ops_cis.o ops_rhf.o ops_io.o ops_mat.o ops_crs.o
OBJ_G20   = grimme-20.o ops_cis.o ops_rhf.o ops_io.o ops_mat.o ops_crs.o

all:    grimme grim20

grimme:	$(OBJ_GRIM)
	$(CC) -o ./bin/grimme $(OBJ_GRIM)  $(LINK_FLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS)

grim20: $(OBJ_G20)
	$(CC) -o ./bin/grimme-20 $(OBJ_G20) $(LINK_FLAGS) $(OPENMP_LIBS) $(LAPACK_LIBS)

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
