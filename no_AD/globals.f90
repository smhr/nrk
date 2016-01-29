module initializeNRK
  implicit none

  integer :: nn != 2051    ! Number of meshpoints
  integer :: ii != 5      ! Number of equations
  integer :: ka != 3      ! Number of boundary conditions at first meshpoint
  integer :: kb != 2       ! Number of boundary conditions at second meshpoint
  double precision :: xb  ! 2nd boundary
  integer :: maxiter != 100 !max num of iterations before giving up
  double precision :: ucy != 1.d0  !convergence speed (do not change unless good reason)
  double precision :: acy != 1.d-7  !desired accuracy of solution
!   integer :: nprevious ! Use previous result (yes=1,no=0)
  
  double precision :: ome2, ome2_up, ome2_low, ome2_step
  double precision :: wave_n, wave_n_up, wave_n_low, wave_n_step
  double precision:: magnetic_field
  
  contains
  
  subroutine readPARAMETER
   implicit none
   open (unit=200, file='nrk.ini', status='old')
   
   read (200,*) nn
   read (200,*) ii
   read (200,*) ka, kb
   read (200,*) xb
   read (200,*) maxiter
   read (200,*) ucy
   read (200,*) acy
!    read (200,*) nprevious
   
   read (200,*) ome2_up, ome2_low, ome2_step
   read (200,*) wave_n_up, wave_n_low, wave_n_step
   read (200,*) magnetic_field
   
  end subroutine readPARAMETER
  
  subroutine printPARAMETER
   implicit none
   open (unit=200, file='nrk.ini', status='old')
   
   write (*,'(a,2x,i10)') 'Number of meshpoints is', nn
   write (*,'(a,2x,i10)') 'Number of equations', ii
   write (*,'(a,2x,2i5)') 'Number of boundary conditions at 1st & 2nd meshpoint is', ka, kb
   write (*,'(a,2x,f8.2)') '2nd boundary is', xb
   write (*,'(a,2x,i5)') 'Max number of iterations before giving result up', maxiter
   write (*,'(a,2x,f6.3)') 'Convergence speed is', ucy
   write (*,'(a,2x,d12.4)') 'Desired accuracy of solution is', acy
!    write (*,'(a,2x,i2)') 'Use previous result (yes=1,no=0)', nprevious
   
   write (*,'(a,2x,3f9.4)') 'Omega^2: Desired max, min & step are', ome2_up, ome2_low, ome2_step
   write (*,'(a,2x,3f9.4)') 'Wave number: Desired max, min & step are:  ', wave_n_up, wave_n_low, wave_n_step
   write (*,'(a,2x,f9.4)') 'magnetic field strength is', magnetic_field
   write (*,*) ''
   
  end subroutine printPARAMETER

end module initializeNRK
