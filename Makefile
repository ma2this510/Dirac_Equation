FC = gfortran
FFLAGS = -O3 -fopenmp -march=native
DEBUG_FLAGS = -g -Wall -Wextra -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=nan -finit-integer=-99999
MARG = Makefile

# MPfun part

mpfuna.o : mpfuna.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfuna.f90

mpfunb.o : mpfunb.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfunb.f90

mpfunc.o : mpfunc.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfunc.f90

mpfund.o : mpfund.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfund.f90

mpfune.o : mpfune.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfune.f90

mpfunf.o : mpfunf.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfunf.f90

mpfung1.o : mpfung1.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfung1.f90

mpfunh1.o : mpfunh1.f90
	$(FC) $(FFLAGS) -ffast-math -c mpfunh1.f90

mpmodule.o : mpmodule.f90
	$(FC) $(FFLAGS) -ffast-math -c mpmodule.f90

mpmask13.o : mpmask13.f90
	$(FC) $(FFLAGS) -ffast-math -c mpmask13.f90

second.o : second.f90
	$(FC) $(FFLAGS) -ffast-math -c second.f90

mpfun : mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung1.o mpfunh1.o mpmodule.o mpmask13.o second.o

# EigenSolver part old

invsg.o : invsg.f90
	$(FC) $(FFLAGS) -c invsg.f90

leq1s.o : leq1s.f90
	$(FC) $(FFLAGS) -c leq1s.f90

# My code part
coefv2.o : coefv2.f90 
	$(FC) $(FFLAGS) -c coefv2.f90

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90

main.out : mpfun invsg.o leq1s.o coefv2.o main.o $(MARG)
	$(FC) $(FFLAGS) -fmax-stack-var-size=0 -o main.out mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung1.o mpfunh1.o mpmodule.o mpmask13.o second.o invsg.o leq1s.o coefv2.o main.o

debug.out : mpfun invsg.o leq1s.o coefv2.f90 main.f90 $(MARG)
	$(FC) $(DEBUG_FLAGS) -c coefv2.f90
	$(FC) $(DEBUG_FLAGS) -c main.f90
	$(FC) $(DEBUG_FLAGS) -o debug.out mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung1.o mpfunh1.o mpmodule.o mpmask13.o second.o invsg.o leq1s.o coefv2.o main.o

run : main.out
	./main.out

clean :
	rm -f *.o *.mod *.out

build : clean main.out

debug : debug.out
