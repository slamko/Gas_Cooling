LIBS=-lraylib -lGL -lm -lpthread -ldl -lrt -lX11 -llapack -lblas

all:
	gfortran -L/home/slamko/.local/lapack -L/home/slamko/src/fortran-raylib -L/home/slamko/proj/gfortran/gas/raylib/lib -Wl,-R/home/slamko/proj/gfortran/gas/raylib/lib -I/home/slamko/src/fortran-raylib gas.f90 /home/slamko/.local/lapack/liblapack.a -lfortran-raylib $(LIBS) -o gas
