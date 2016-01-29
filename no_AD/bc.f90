!***********************************************************
subroutine initializev(v)
  use initializeNRK
  implicit none

  integer, dimension(ii) :: v
  integer :: i
  
! *** Enter the permutation v-vector for the first boundary condition 
! *** matrix. Ask me for detail.  

  do i=1,ii
     v(i) = i
  enddo
  
  v(2) = 4
  v(4) = 2
  
end subroutine initializev


! ***********************************************************************

subroutine bc(xxa,xxb,ya,yb,g,gs)
  use initializeNRK
  implicit none
  
  double precision :: xxa,xxb 
  double precision, dimension(ii) :: ya, yb, g
  double precision, dimension(ii,ii) :: gs
  integer :: i
  
! ********************************************************************
! *** Subroutine in which the boudary condition functions g must
!     be entered. The only constraint is that the bcs at the first 
!     boundary must be entered first, and those at the second boundary 
!     must be entered afterwards.
! *** Input parameters:
!     xa, xb: first and last meshpoints 
!     ya,yb: vectors of guesses at the first and last meshpoints
!     g: functions g(x,y) evaluated at meshpoints considered.
!     gs: Jacobian of g

! Example presented here has: 
! y(xa) = 0
! dy/dx(xa) = 1 
! y(xb) = 0
! 
!*********************************************************************

!      write(*,*) '- Starting bc - '

!begin_smhr
! !   g(1) = ya(1)
! !   gs(1,1) = 1.d0
! !   
! !   g(2) = ya(2) - 1.d0
! !   gs(2,2) = 1.d0
! ! 
! !   g(3) = yb(1)
! !   gs(3,1) = 1.d0

!   do i = 1,ii
!     write(6,1) i, ya(i), yb(i), ya(i), yb(i)
!   1 format ('i,ya(i),yb(i)',i3,2e10.2,2f10.3)
!   enddo
  
  g(1) = ya(1) - 1.
  gs(1,1) = 1.d0
  
  g(2) = ya(4)
  gs(2,4) = 1.d0
  
  g(3) = ya(3)
  gs(3,3) = 1.d0

  g(4) = yb(1)
  gs(4,1) = 1
  
  g(5) = yb(3)
  gs(5,3) = 1
  
!   do i = 1,ii
!     write(6,3) i,g(i),g(i)
!   3 format ('i,g(i)',i3,e10.2,f8.3)
!   enddo
  
!end_smhr

!      write(*,*) '- Starting bc - '
  
end subroutine bc

! *****************************************************************





