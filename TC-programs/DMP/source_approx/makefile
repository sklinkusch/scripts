.SUFFIXES:
.SUFFIXES: .f .f90 .o


# if one compiles on Itanium II
FL = gfortran 
CMPLFLG = -c -O3 
LNKFLG = -o $(TARGET) -llapack -lblas 
OPENMPFLG = -fopenmp
#CMPLFLG = -c -march=i386 -funroll-loops -O3 -msse3 -ftree-vectorize
#LNKFLG = -o $(TARGET) /Users/jct/bin/lapack-3.2.1/lapack_OSX.a /Users/jct/bin/GotoBLAS/libgoto.a

TARGET = rhoTDCI.x

OBJS = scalarOP.o matrix.o rungekutta.o spline.o rhoprop.o

$(TARGET) : $(OBJS)
	$(FL) $(OBJS) $(LNKFLG) $(OPENMPFLG)

#%.f : %.o
.f90.o:
	$(FL) $(CMPLFLG) $(OPENMPFLG) $<

.f.o:
	$(FL) $(CMPLFLG) $(OPENMPFLG) $<

#
#  clean	- remove all object files
#
clean:
	rm -f *.o *.mod  core
