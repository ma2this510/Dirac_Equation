FC = gfortran
FFLAGS = -llapack -lblas -O3 -fopenmp -ffast-math -march=native -fdefault-real-8
DEBUG_FLAGS = -llapack -lblas -Wall -Og -g -fdefault-real-8

MARG = Makefile

coefv2.o : coefv2.f90 
	$(FC) $(FFLAGS) -c coefv2.f90

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90

main.out : coefv2.o main.o $(MARG)
	$(FC) $(FFLAGS) -o main.out coefv2.o main.o

debug.out : coefv2.f90 main.f90 $(MARG)
	$(FC) $(DEBUG_FLAGS) -c coefv2.f90
	$(FC) $(DEBUG_FLAGS) -c main.f90
	$(FC) $(DEBUG_FLAGS) -o debug.out coefv2.o main.o

run : main.out
	./main.out

clean :
	rm -f *.o *.mod *.out

build : clean main.out

debug : clean debug.out
