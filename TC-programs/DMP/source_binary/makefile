.SUFFIXES:
.SUFFIXES: .f .f90 .o


# if one compiles on Itanium II
FL = gfortran 
CMPLFLG = -c -O3 -Wall
LNKFLG = -o $(TARGET) -O3 -Wall -llapack -lblas
OPENMPFLG = -fopenmp
#CMPLFLG = -c -O2 -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics
#LNKFLG = -o $(TARGET) -llapack -lblas -O2 -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -fall-intrinsics

TARGET = rhoTDCI.bin

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
