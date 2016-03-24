mkl77 = disable
#mkl77 = enable
loop_mod = off
brf_mod = on

FC = ifort
# FC = gfortran
ifeq ($(FC), gfortran)
 FFLAGS = -fbounds-check -fbacktrace
else 
 FFLAGS = 
endif
ifeq ($(mkl77), enable)
 loop_mod = off
 brf_mod = on
 MKLROOT= /media/d3/Intel/composer_xe_2013.3.163/mkl
 MKLOPTIONS = -i8 -I$(MKLROOT)/include -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl
 FFLAGS += $(MKLOPTIONS)
 BINNAME = NRK_AD_brf_no_loop_mkl.exe
 SOURCE = globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 lapacky.f90 rhs.f90 ODEsolve_brf_no_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90
endif

ifeq ($(mkl77), disable)
 ifeq ($(loop_mod), on)
  ifeq ($(brf_mod), on)
   BINNAME = NRK_AD_brf_loop.exe
   SOURCE = globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_brf_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90
  endif
 else
  ifeq ($(brf_mod), on)
   BINNAME = NRK_AD_brf_no_loop.exe
   SOURCE = globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_brf_no_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90
  else
   BINNAME = NRK_AD_no_loop.exe
   SOURCE = globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_no_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90
  endif
 endif
endif
  
NRK:
	$(FC) $(SOURCE) $(FFLAGS) -o $(BINNAME)
# NRK_AD_brf_loop:
# 	gfortran -g -fbounds-check -fbacktrace globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_brf_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90 -o NRK_AD_brf_loop.exe
# NRK_AD_brf_no_loop:
# # 	gfortran globals.f90 bc.f90 guess.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve.f90 amean.f90 lhs.f90 print.f90 -o NRK.exe
# 	gfortran -g -fbounds-check -fbacktrace globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_brf_no_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90 -o NRK_AD_brf_no_loop.exe
# 
# NRK_AD_no_loop:
# 	gfortran -g -fbounds-check -fbacktrace globals.f90 bc.f90 density.f90 guess.f90 main.f90 mesh.f90 mynrk.f90 rhs.f90 ODEsolve_no_loop.f90 amean.f90 lhs.f90 print.f90 sweep_blanks.f90 -o NRK_AD_no_loop.exe

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
	rm -f ./*.o ./*.mod ./*.exe ./*.dat ./final_res 

