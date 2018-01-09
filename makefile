CC = icc -openmp -xSSE4.2 -O3 -ansi-alias -qopt-report=1
FC = ifort -openmp -xSSE4.2 -O3 -ansi-alias  -qopt-report=1
INC=
CFLAGS= -std=c99 -I/opt/intel/mkl/include/ -I./src/
LFLAGS= -L/opt/intel/mkl/lib/intel64/ -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lgfortran -lifcore

all:
	cd src && make && cd ..
gaussian: all
	rm -rf gaussian
	rm -rf gaussian.o
	${CC} ${CFLAGS} -I./src -c gaussian.c
	${CC} ${LFLAGS} -o gaussian gaussian.o ./src/libMKJ.a
example: all
	rm -rf example
	rm -rf example.o
	${CC} ${CFLAGS} -I./src -c example.c
	${CC} ${LFLAGS} -o example example.o ./src/libMKJ.a
metal: all
	rm -rf metal
	rm -rf metal.o
	${CC} ${CFLAGS} -I./src -c metal.c
	${CC} ${LFLAGS} -o metal metal.o ./src/libMKJ.a
Fig6: all
	rm -rf Fig6
	rm -rf Fig6.o
	${CC} ${CFLAGS} -I./src -c Fig6.c
	${CC} ${LFLAGS} -o Fig6 Fig6.o ./src/libMKJ.a
Fig9: all
	rm -rf Fig9
	rm -rf Fig9.o
	${CC} ${CFLAGS} -I./src -c Fig9.c
	${CC} ${LFLAGS} -o Fig9 Fig9.o ./src/libMKJ.a
clean:
	cd src && make clean && cd ..
	rm -rf Fig9.o
	rm -rf Fig9