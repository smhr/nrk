subroutine density (rho_val)
use initializeNRK
  implicit none
  
  double precision :: xa = 0. ! Position of first boundary
!   double precision :: xb = 10.d0 ! Postition of second boundary

  double precision, dimension(nn) :: x
  double precision, dimension(nn),intent(out) :: rho_val
  integer :: i

! *** ***************************************************************
! *** User-supplied subroutine in which the mesh is entered
! *** Input parameters:
!     x : nn-vector of meshpoints (needs to be strictly monotonic, but
!     can either increase or decrease).
!     nn: total number of meshpoints.
! ***************************************************************

!   write(*,*) '- Starting mesh -'
  open (unit=91, file='rho.dat')
  do i=1,nn
     x(i) = xa+(xb-xa)*(i-1.d0)/(nn-1.d0)
     rho_val(i) = 1.d0/((1. + x(i)*x(i)/8.d0)**2.d0)
     write(91,*) x(i), rho_val(i)
  enddo
!   flush(91)
  
!   write(*,*) '- Ending mesh -'
  
end subroutine density