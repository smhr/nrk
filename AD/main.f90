program main
 use initializeNRK
 implicit none
  call readPARAMETER ()
  call ODEsolve ()
end program main