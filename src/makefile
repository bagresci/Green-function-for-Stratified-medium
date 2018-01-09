CC = icc -openmp -xSSE4.2 -O3 -ansi-alias -qopt-report=5
FC = ifort -openmp -xSSE4.2 -O3 -ansi-alias -qopt-report=5
CFLAGS= -std=c99 -I/opt/intel/mkl/include/ -I./src/
LFLAGS= -lhdf5 -L/opt/intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_sequential
TARGET = libMKJ.a
OBJECT = gauss_legendre.o bfor.o zbesh.o bessel.o machine.o func.o mem.o scattering.o

all: ${OBJECT}
	ar rcs ${TARGET} ${OBJECT}

gauss_legndre.o: gauss_legendre.c gauss_legendre.h
	${CC} ${CFLAGS} -c gauss_legendre.c
bfor.o: bfor.f90 machine.for zbesh.for
	${FC} -c bfor.f90
zbesh.o: zbesh.for machine.for
	${FC} -c zbesh.for
machine.o: machine.for
	${FC} -c machine.for
bessel.o: bessel.c bfor.o
	${CC} ${CFLAGS} -c bessel.c
func.o: func.c mkj.h
	${CC} ${CFLAGS} -c func.c
mem.o: mem.c mkj.h
	${CC} ${CFLAGS} -c mem.c
scattring.o: scattering.c func.o mkj.h 
	${CC} ${CFLAGS} -c scattering.c

clean:
	rm -rf ${OBJECT} ${TARGET}
