subroutine mesh(x,nprevious)
use initializeNRK
  implicit none
  
  double precision :: xa = 0.d0 ! Position of first boundary
!   double precision :: xb = 10.d0 ! Postition of second boundary

  double precision, dimension(nn),intent(out) :: x
!   double precision, dimension(nn),intent(out) :: rho_val
  integer :: i, nprevious

! *** ***************************************************************
! *** User-supplied subroutine in which the mesh is entered
! *** Input parameters:
!     x : nn-vector of meshpoints (needs to be strictly monotonic, but
!     can either increase or decrease).
!     nn: total number of meshpoints.
! ***************************************************************

!   write(*,*) '- Starting mesh -'
if(nprevious.eq.2) then
  
    open(12,file='result.dat.org',status='old')
     
    do i=1,nn
       read(12,*) x(i)
    enddo
    close(12)
elseif (nprevious.eq.0) then
  do i=1,nn
     x(i) = xa+(xb-xa)*(i-1.d0)/(nn-1.d0)
!      rho_val(i) = 1./((1. + x(i)*x(i)/8.)**2.)
!      write(*,*) x(i)
  enddo
endif
end subroutine mesh
