FC = gfortran
FFLAGS = -O3 -fopenmp -march=native -ffast-math
DEBUG_FLAGS = -fopenmp -g -Wall -Wextra -ffpe-trap=invalid,zero,overflow -fbounds-check -finit-real=nan -finit-integer=-99999
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

# EigenSolver part new

pythag.o : pythag.f90
	${FC} ${FFLAGS} -c pythag.f90

rebak.o : rebak.f90
	${FC} ${FFLAGS} -c rebak.f90

reduc.o : reduc.f90
	${FC} ${FFLAGS} -c reduc.f90

rsg.o : rsg.f90
	${FC} ${FFLAGS} -c rsg.f90

tql2.o : tql2.f90
	${FC} ${FFLAGS} -c tql2.f90

tqlrat.o : tqlrat.f90
	${FC} ${FFLAGS} -c tqlrat.f90

tred1.o : tred1.f90
	${FC} ${FFLAGS} -c tred1.f90

tred2.o : tred2.f90
	${FC} ${FFLAGS} -c tred2.f90

eigen : pythag.o rebak.o reduc.o rsg.o tql2.o tqlrat.o tred1.o tred2.o

# My code part
tools_mp.o : tools_mp.f90
	$(FC) $(FFLAGS) -c tools_mp.f90

bspline_gen.o : bspline_gen.f90
	$(FC) $(FFLAGS) -c bspline_gen.f90

coefv2.o : coefv2.f90 
	$(FC) $(FFLAGS) -c coefv2.f90

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90

main.out : mpfun eigen tools_mp.o bspline_gen.o coefv2.o main.o $(MARG)
	$(FC) $(FFLAGS) -fmax-stack-var-size=0 -o main.out mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung1.o mpfunh1.o mpmodule.o mpmask13.o second.o pythag.o rebak.o reduc.o rsg.o tql2.o tqlrat.o tred1.o tred2.o tools_mp.o bspline_gen.o coefv2.o main.o

debug.out : mpfun eigen tools_mp.o bspline_gen.o coefv2.f90 main.f90 $(MARG)
	$(FC) $(DEBUG_FLAGS) -c coefv2.f90
	$(FC) $(DEBUG_FLAGS) -c main.f90
	$(FC) $(DEBUG_FLAGS) -o debug.out mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung1.o mpfunh1.o mpmodule.o mpmask13.o second.o pythag.o rebak.o reduc.o rsg.o tql2.o tqlrat.o tred1.o tred2.o tools_mp.o bspline_gen.o coefv2.o main.o

run : main.out
	./main.out

clean :
	rm -f *.o *.mod *.out

build : clean main.out

debug : debug.out
