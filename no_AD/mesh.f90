subroutine mesh(x,mesh_selector)
use initializeNRK
  implicit none
  
  double precision :: xa = 0. ! Position of first boundary
!   double precision :: xb = 10.d0 ! Postition of second boundary

  double precision, dimension(nn),intent(out) :: x
!   double precision, dimension(nn),intent(out) :: rho_val
  integer :: i

! *** ***************************************************************
! *** User-supplied subroutine in which the mesh is entered
! *** Input parameters:
!     x : nn-vector of meshpoints (needs to be strictly monotonic, but
!     can either increase or decrease).
!     nn: total number of meshpoints.
! ***************************************************************

!   write(*,*) '- Starting mesh -'
  if (mesh_selector == 0) then
     do i=1,nn
        x(i) = xa+(xb-xa)*(i-1.d0)/(nn-1.d0)
     enddo
  elseif (mesh_selector == 1)
     open (10, file='mesh.dat')
     do i=1,nn
        read (10,*) x(i)
     enddo
  endif
!   flush(91)
  
!   write(*,*) '- Ending mesh -'
  
end subroutine mesh
