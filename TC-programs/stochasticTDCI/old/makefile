.SUFFIXES:
.SUFFIXES: .f .f90 .o

FL = gfortran
CMPLFLG = -c -O3 -fopenmp
#CMPLFLG = -c -g -fcheck=all -Wall -fbacktrace  ### extremely useful debugging tools
LNKFLG = -o $(TARGET) -llapack -lblas -fopenmp #-lgomp

TARGET = sTDCIold.x

OBJS = localfield.o matrix.o rungekutta.o \
       stepper.o stoch_prop.o

$(TARGET) : $(OBJS)
	$(FL) $(OBJS) $(LNKFLG)

#%.f : %.o
.f90.o:
	$(FL) $(CMPLFLG) $<

.f.o:
	$(FL) $(CMPLFLG) $<

#
#  clean	- remove all object files
#
clean:
	rm -f *.o *.mod  core
