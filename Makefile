compiler = gfortran
debug = false
profiling = false
linux = false

ifeq ($(debug),true)
	flag1 = -O0 -fbacktrace -pedantic -fall-intrinsics -Wconversion -pedantic
	flag2 = -Wall -Wextra -fopenmp -fcheck=all -ffree-form -finit-real=nan -fdump-core
	flag3 = -ffpe-trap=invalid,zero,overflow,underflow  -fimplicit-none -fbounds-check
	flag = $(flag1) $(flag2) $(flag3)
else
	flag = -O3 -fopenmp
endif

ifeq ($(profiling),true)
	flags = $(flag) -p
else
	flags = $(flag)
endif

ifeq ($(linux),true)
	lib = lib/libcfitsio.so.3
else
	lib = lib/libcfitsio.5.dylib
endif

all: init ARTES clean

init:
	@echo \####################################################
	@echo \ \ \ \ \ \ \ \ \ \ \  _ \   \  ___ \  _____ \ \ ___ \ \ ___ 
	@echo \ \ \ \ \ \ \ \ \ \   \/_\\  \  \| _ \\ \|_ \  _\| \| __\| \/ __\|
	@echo \ \ \ \ \ \ \ \ \  \/ _ \\ \ \| \  \/ \  \| \| \  \| _\| \ \\__ \\
	@echo \ \ \ \ \ \ \ \ \ \/_\/ \\_\\ \|_\|_\\ \  \|_\| \  \|___\| \|___\/
	@echo
	@echo Atmospheric Radiative Transfer for Exoplanet Science
	@echo
	@echo \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ Developed\ by:
	@echo
	@echo \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ Tomas Stolker
	@echo \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ T.Stolker\ [at]\ uva.nl
	@echo
	@echo \ \ \ \ \ \ Anton\ Pannekoek\ Institute\ for\ Astronomy
	@echo
	@echo ----------------------------------------------------
	@echo
	@echo \-\-\>\ Initialize
	@echo
	mkdir -p input/
	mkdir -p output/
	mkdir -p bin/
	@echo
	@echo ----------------------------------------------------
	@echo

ARTES: $(ARTES)
	@echo \-\-\>\ Compile
	@echo
	$(compiler) $(flags) -c src/ARTES.f90
	@echo
	$(compiler) $(flags) -o bin/ARTES ARTES.o $(lib)
	@echo
	@echo ----------------------------------------------------
	@echo

clean:
	@echo \-\-\>\ Cleanup
	@echo
	rm -f *.o
	rm -f *.mod
	rm -f gmon.out
	@echo
	@echo \####################################################
