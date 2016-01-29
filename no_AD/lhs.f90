! ***********************************************************************
! *** Subroutines where the matrix M_ij and its Jacobian M_ijl are entered. 
!     The ones here are the default; Only modify when M are non-unity.
! ***********************************************************************

double precision function am(i,j,x,y,in)
  use initializeNRK
  implicit none
  
  ! Input: indices i and j
  ! Input: mesh x (by definition, x at the meshpoint considered)
  ! Input: function array at meshpoint considered
  ! Input: in = index of meshpoint considered 

  double precision, dimension(ii) :: y
  double precision :: x, rho, pd, pdd, b
  integer :: i,j,in
!begin_smhr 
  pd = 1. ; pdd = 0.
  b = magnetic_field
  rho = 1./((1. + x*x/8.)**2.)
!   print*,'x,rho',x,rho
  am=0.d0
  if (i==j) am = 1.
  if (i==1.and.j==1) am = x*(pd+b*b/rho*(1.-y(5)*y(5)/ome2*pd))
  if (i==1.and.j==2) am = x*(rho-y(5)*y(5)*b*b/ome2)
  if (i==3.and.j==3) am = x*rho
  if (i==4.and.j==4) am = x
  if (i==4.and.j==2) am = 1.
  
!   write (*,*) "##", i, j, am
  
!end_smhr
  
  
!   if(i.eq.j) am=1.

end function am

! **********************************************************************

double precision function amd(i,j,l,x,y,in)
  use initializeNRK
  implicit none

  ! Input: indices i and j
  ! Input: mesh x (by definition, x at the meshpoint considered)
  ! Input: function array at meshpoint considered
  ! Input: in = index of meshpoint considered 
  
  double precision, dimension(ii) :: y
  double precision :: x, b, pd, pdd, rho
  integer :: i,j,l,in
  
  pd = 1. ; pdd = 0.
  b = magnetic_field
  rho = 1./((1. + x*x/8.)**2.)
  
  amd=0.d0
  if (i==1.and.j==1.and.l==5) amd = -2*x*b*b/(ome2*rho)*pd*y(5)
  if (i==1.and.j==2.and.l==5) amd = -2*x*b*b/(ome2)*y(5)
  
end function amd

! ************************


