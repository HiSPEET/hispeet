# This is a Makefile, recipes are annoted here and made with make

.phony: all, clean, tidy

PROGRAM_NAMES= ping_pong sem_1d
PROGRAMS     = $(addprefix bin/,$(PROGRAM_NAMES))

all:    $(PROGRAMS)

clean:
	rm -f *.{mod,o} tmp/*.{mod,o}

tidy:   clean
	rm -f $(PROGRAMS)

# Compiler & flags -------------------------------------------------------------

FC     := mpifort
FFLAGS := -O3 -ffree-form -fopenmp -Wall -std=gnu -J tmp/ -I tmp/

# Rules ------------------------------------------------------------------------

# generate binary, use all available object files
bin/%: tmp/%.o tmp/kind_parameters.o tmp/constants.o
	$(FC) $(FFLAGS) $^ -o $@

# generate object files from source files
tmp/%.o: src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Module dependencies ----------------------------------------------------------

tmp/constants.o:                                                               \
        tmp/kind_parameters.o

tmp/gauss_jacobi.o:                                                            \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o

tmp/ping_pong.o:                                                               \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o

tmp/sem_1d.o:                                                                  \
        tmp/kind_parameters.o                                                  \
        tmp/constants.o                                                        \
        tmp/gauss_jacobi.o

# Program dependencies .........................................................

bin/ping_pong:

bin/sem_1d:                                                                    \
	tmp/gauss_jacobi.o
