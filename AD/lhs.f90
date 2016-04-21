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
  double precision :: x, rho, pd, pdd, b, j1, eh
  integer :: i,j,in
!begin_smhr 
  eh = eta
  if (x > ADlimit) eh = 0.d0
  pd = 1.d0 ; pdd = 0.d0
  rho = 1.d0/((1.d0 + x*x/8.d0)**2.d0)
  b = magnetic_field
!  b = magnetic_field * rho

  j1 = -b*b*b*y(5)*y(5)*eh/(eh*b*b*y(5)*y(5) + omegg)
!   print*,'x,rho',x,rho
  am=0.d0
  if (i==j) am = 1.d0
  if (i==1.and.j==1) am = pd
  if (i==1.and.j==2) am = rho
  if (i==1.and.j==6) am = j1+b
  if (i==3.and.j==3) am = x*rho
  if (i==4.and.j==4) am = x
  if (i==4.and.j==2) am = 1.d0
  if (i==7.and.j==6) am = eh*b*b
  if (i==7.and.j==7) am = eh*b*b*x
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
  double precision :: x, b, pd, pdd, rho, eh
  integer :: i,j,l,in
  
  eh = eta
  if (x > ADlimit) eh = 0.d0
  pd = 1.d0 ; pdd = 0.d0
  rho = 1.d0/((1.d0 + x*x/8.d0)**2.d0)
  b = magnetic_field
!  b = magnetic_field * rho
  
  amd=0.d0
  if (i==1.and.j==6.and.l==5) &
   & amd = (-2.d0*b**3.d0*eh*y(5)*(eh*b*b*y(5)*y(5) + omegg) &
         & -2.d0*b*b*eh*y(5)*(-b*b*b*eh*y(5)*y(5))) &
         & /(eh*b*b*y(5)*y(5) + omegg)**2.d0
  
end function amd

! ************************


