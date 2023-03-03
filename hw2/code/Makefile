FC = gfortran
FFLAGS = -O3 -Wall -fcheck=all -fbacktrace -g -finit-real=NaN
SRC = mod_tools.f90 mod_hw2.f90 main.f90
SRC_EXTRA = mod_tools.f90 mod_hw2.f90 extra.f90
OBJ = ${SRC:.f90=.o}
OBJ_EXTRA = ${SRC_EXTRA:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

extra: $(OBJ_EXTRA)
	$(FC) $(FFLAGS) -o $@ $(OBJ_EXTRA)

clean:
	@rm -f *.mod *.o main extra