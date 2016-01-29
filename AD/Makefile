NRK:
# 	gfortran globals.f90 bc.f90 guess.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve.f90 amean.f90 lhs.f90 print.f90 -o NRK.exe
	gfortran -g -fbounds-check -fbacktrace globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90 -o NRK.exe

# FFLAGS = -fbounds-check -fbacktrace
# FC = gfortran
# SOURCE = globals.f90 bc.f90 guess.f90 mesh.f90 mynrk.f90 \
# rhs.f90 ODEsolve.f90 amean.f90 lhs.f90 print.f90
# OBJECTS = $(SOURCE:.f90=.o)
# 
# NRK: initializenrk.mod $(OBJECTS)
# 	$(FC) $(OBJECTS) $(FFLAGS) -o NRK.exe
# initializenrk.mod: globals.f90
# 	$(FC) $(FFLAGS) -c globals.f90
clean:
	rm -f ./*.o ./*.mod ./NRK.exe ./*.dat ./final_res 
