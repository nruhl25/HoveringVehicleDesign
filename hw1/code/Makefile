FC = gfortran
FFLAGS = -O3 -Wall -fcheck=all -fbacktrace -g -finit-real=NaN
SRC = mod_copter.f90 mod_splines.f90 mod_tools.f90 main.f90
OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o main