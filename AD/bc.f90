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
  
!   v(2) = 4
!   v(4) = 2

!   v(2) = 4
!   v(4) = 7
!   v(7) = 2
 if (bcType == 431) then
    v(2) = 4
    v(4) = 6
    v(6) = 2
 elseif (bcType == 341) then
    v(2) = 4
    v(4) = 2
! elseif (bcType == 345) then
!   v(2) = 3
!   v(3) = 4
!   v(4) = 2
 endif
  
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (bcType == 431) then
    g(1) = ya(1) - 1.d0
    gs(1,1) = 1.d0
  
    g(2) = ya(4)
    gs(2,4) = 1.d0
  
    g(3) = ya(3)
    gs(3,3) = 1.d0
  
!   g(4) = ya(7) - (-0.1573847)
!   gs(4,7) = 1.
  
!   g(4) = ya(6) - (0.0279609)
    g(4) = 0.d0
    gs(4,6) = 1.d0
  
    g(5) = yb(1)
    gs(5,1) = 1.d0
  
    g(6) = yb(2)
    gs(6,2) = 1.d0
  
    g(7) = yb(3)
    gs(7,3) = 1.d0
  
  elseif (bcType == 432) then

    g(1) = ya(1) - 1.d0
    gs(1,1) = 1.d0
  
    g(2) = ya(4)
    gs(2,4) = 1.d0
  
    g(3) = ya(3)
    gs(3,3) = 1.d0
  
!   g(4) = ya(7) - (-0.1573847)
!   gs(4,7) = 1.
  
!   g(4) = ya(6) - (0.0279609)
    g(4) = 0.d0
    gs(4,6) = 1.d0
  
    g(5) = yb(1)
    gs(5,1) = 1.d0
  
    g(6) = yb(2)
    gs(6,2) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  
  elseif (bcType == 341) then
     g(1) = ya(1) - 1.d0
     gs(1,1) = 1.d0
  
     g(2) = ya(4)
     gs(2,4) = 1.d0
  
     g(3) = ya(3)
     gs(3,3) = 1.d0
  
     g(4) = yb(6)
     gs(4,6) = 1.d0
  
     g(5) = yb(1)
     gs(5,1) = 1.d0
  
     g(6) = yb(2)
     gs(6,2) = 1.d0
  
     g(7) = yb(3)
     gs(7,3) = 1.d0
  elseif (bcType == 342) then
     g(1) = ya(1) - 1.d0
     gs(1,1) = 1.d0
  
     g(2) = ya(4)
     gs(2,4) = 1.d0
  
     g(3) = ya(3)
     gs(3,3) = 1.d0
  
     g(4) = yb(6)
     gs(4,6) = 1.d0
  
     g(5) = yb(1)
     gs(5,1) = 1.d0
  
     g(6) = yb(7)
     gs(6,7) = 1.d0
  
     g(7) = yb(3)
     gs(7,3) = 1.d0
  elseif (bcType == 343) then
     g(1) = ya(1) - 1.d0
     gs(1,1) = 1.d0
  
     g(2) = ya(4)
     gs(2,4) = 1.d0
  
     g(3) = ya(3)
     gs(3,3) = 1.d0
  
     g(4) = 0.d0
     gs(4,6) = 1.d0
  
     g(5) = yb(1)
     gs(5,1) = 1.d0
  
     g(6) = yb(6)
     gs(6,6) = 1.d0
  
     g(7) = yb(3)
     gs(7,3) = 1.d0

  elseif (bcType == 344) then

    g(1) = ya(1) - 1.d0
    gs(1,1) = 1.d0
  
    g(2) = ya(4)
    gs(2,4) = 1.d0
  
    g(3) = ya(3)
    gs(3,3) = 1.d0
  
!   g(4) = ya(7) - (-0.1573847)
!   gs(4,7) = 1.
  
!   g(4) = ya(6) - (0.0279609)
    g(4) = yb(3)
    gs(4,3) = 1.d0
  
    g(5) = yb(1)
    gs(5,1) = 1.d0
  
    g(6) = yb(2)
    gs(6,2) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  
  elseif (bcType == 345) then

    g(1) = ya(1) - 1.d0
    gs(1,1) = 1.d0
  
    g(2) = ya(3)
    gs(2,3) = 1.d0
  
    g(3) = ya(4)
    gs(3,4) = 1.d0
  
    g(4) = yb(1)
    gs(4,1) = 1.d0
  
    g(5) = yb(2)
    gs(5,2) = 1.d0
  
    g(6) = yb(6)
    gs(6,6) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  
  elseif (bcType == 161) then

    g(1) = ya(3)
    gs(1,3) = 1.d0
  
    g(2) = yb(1)
    gs(2,1) = 1.d0
  
    g(3) = yb(2)
    gs(3,2) = 1.d0
  
    g(4) = yb(3)
    gs(4,3) = 1.d0
  
    g(5) = yb(4)
    gs(5,4) = 1.d0
  
    g(6) = yb(6)
    gs(6,6) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  
  elseif (bcType == 251) then

    g(1) = ya(3)
    gs(1,3) = 1.d0
  
    g(2) = ya(4)
    gs(2,4) = 1.d0
  
    g(3) = yb(1)
    gs(3,1) = 1.d0
  
    g(4) = yb(3)
    gs(4,3) = 1.d0
  
    g(5) = yb(2)
    gs(5,2) = 1.d0
  
    g(6) = yb(6)
    gs(6,6) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  elseif (bcType == 346) then

    g(1) = ya(3)
    gs(1,3) = 1.d0
  
    g(2) = ya(4)
    gs(2,4) = 1.d0
  
    g(3) = 0.d0
    gs(3,5) = 1.d0

    g(4) = yb(3)
    gs(4,3) = 1.d0
  
    g(5) = 0.d0
    gs(5,5) = 1.d0
  
    g(6) = yb(6)
    gs(6,6) = 1.d0
  
    g(7) = yb(7)
    gs(7,7) = 1.d0
  
 elseif (bcType == 347) then

    g(1) = ya(1) - 1.d0
    gs(1,1) = 1.d0
  
    g(2) = ya(3)
    gs(2,3) = 1.d0
  
    g(3) = ya(4)
    gs(3,4) = 1.d0
  
    g(4) = yb(1)
    gs(4,1) = 1.d0
  
    g(5) = yb(2)
    gs(5,2) = 1.d0
  
    g(6) = yb(3)
    gs(6,3) = 1.d0
  
    g(7) = yb(6)
    gs(7,6) = 1.d0
  
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   g(1) = ya(1) - 1.
!   gs(1,1) = 1.d0
!   
!   g(2) = ya(4)
!   gs(2,4) = 1.d0
!   
!   g(3) = ya(3)
!   gs(3,3) = 1.d0
!   
!   g(4) = ya(2) + ya(1)/(ya(5)*ya(5))
!   gs(4,2) = 1.
!   gs(4,1) = 1./(ya(5)*ya(5))
!   gs(4,5) = - 2.*ya(1)/(ya(5)**3.)
! 
!   g(5) = yb(1)
!   gs(5,1) = 1.
!   
!   g(6) = yb(2)
!   gs(6,2) = 1.
!   
!   g(7) = yb(3)
!   gs(7,3) = 1.
  
!   g(7) = yb(4)
!   gs(7,4) = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!  
!   g(1) = ya(1) - 1.
!   gs(1,1) = 1.d0
!   
!   g(2) = ya(4)
!   gs(2,4) = 1.d0
!   
!   g(3) = ya(3)
!   gs(3,3) = 1.d0
!   
! !   g(4) = ya(2) + ya(1)/(ya(5)*ya(5))
! !   gs(4,2) = 1.
! !   gs(4,1) = 1./(ya(5)*ya(5))
! !   gs(4,5) = - 2.*ya(1)/(ya(5)**3.)
! 
!   g(4) = yb(1)
!   gs(4,1) = 1.
!   
!   g(5) = yb(2)
!   gs(5,2) = 1.
!   
!   g(6) = yb(3)
!   gs(6,3) = 1.
!   
!   g(7) = yb(4)
!   gs(7,4) = 1.
!!!!!!!!!!!!!!!!!!!!!!!!!  
!   do i = 1,ii
!     write(6,3) i,g(i),g(i)
!   3 format ('i,g(i)',i3,e10.2,f8.3)
!   enddo
  
!end_smhr

!      write(*,*) '- Starting bc - '
  
end subroutine bc

! *****************************************************************





