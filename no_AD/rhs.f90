! **********************************************************************
subroutine rhs(x,y,f,fd,i)
  use initializeNRK
  implicit none
      
  double precision, dimension(ii) :: y, f
  double precision, dimension(ii,ii) :: fd
  double precision :: rho, rhod, rhodd, psid, A1, A2, A3, pd, pdd, b, c1, k, c2
  double precision :: x
  integer :: i

! *** Subroutine in which the RHS of the system is entered, as
!     well as its Jacobian. User provided.
! 
! Input parameters:
!     x = point considered
!     Y = vector of functions at point considered
!     f(k) = f_k at point considered
!     fd(k,j) = df_k/dy_j at point considered 
!     i = number of current meshpoint (in main program,x(i) = x)
!     
! ******************************************************************
      
!      write(*,*) '- Starting rhs - '

! To be completed by user: F_i function (f) and Jacobian (fd).
! if Jacobian matrix component is zero - no need to enter it 

! Here is an example on y'' = -omega**2 y as in class.
! y(1) = y
! y(2) = dy/dx
! y(3) = omega^2  (the eigenvalue)

!begin_smhr
!   f(1) = y(2)
!   fd(1,2) = 1.d0
!   
!   f(2) = -y(3)*y(1)
!   fd(2,1) = -y(3)
!   fd(2,3) = -y(1)
!   
!   f(3) = 0.d0
!   k = 0.5
!   ome2 = -0.04
!   print*,'x=',x

  pd = 1. ; pdd = 0.
  b = magnetic_field
  k = y(5)
  rho = 1./((1. + x*x/8.)**2.)
  rhod = -2. * (x/4.) / ((1. + x*x/8.)**3.)
  rhodd = 3./8.*x*x / ((1. + x*x/8.)**4.) - 0.5/((1. + x*x/8.)**3.)
  psid = x/(2.+x*x/4.)

  ! matrix f coefficients
!   c1 = (k*k*b*b)/ome2
  c1 = (k*k*b*b)
!   a1 = x*(-c1/rho*(pdd-2.*pd/rho)*rhod + (pdd-2.*b*b/(rho*rho))*rhod + psid)
  a1 = x*(-c1/rho*(pdd-2.*pd/rho)*rhod + ome2*(pdd-2.*b*b/(rho*rho))*rhod + ome2*psid)
  a2 = x*(c1/rho * rhod)
!   a3 = x*b*b/ome2 * (k*k + rhodd/rho - 2.*(rhod/rho)**2.) - b*b/ome2 * rhod/rho - x*rho
  a3 = x*b*b * (k*k + rhodd/rho - 2.*(rhod/rho)**2.) - b*b * rhod/rho - ome2*x*rho
  
  ! matrix fd coefficients
!   c2 = (2.*k*b*b)/ome2
  c2 = (2.*k*b*b)
  
!   write(6,'(a25,i8,4f8.3)')'i,x,rho,rhod,psid',i,x,rho,rhod,psid
    
!   f(1) = y(3)*rho - y(1)*psid
  f(1) = -a1*y(1) - a2*y(2) - a3*y(3)
!   fd(1,1) = -psid
  fd(1,1) = -a1
  fd(1,2) = -a2
!   fd(1,3) = rho
  fd(1,3) = -a3
  fd(1,5) = x*(c2/rho*(pdd-2.*pd/rho)*rhod * y(1) - (c2/rho * rhod) * y(2) - c2 * y(3))
    
  f(2) = y(4)
  fd(2,4) = 1.
   
  f(3) = -y(3)*rho + x*(y(5)*y(5)-ome2)*y(1) + x*y(5)*y(5)*rho*y(2) - x*y(3)*rhod
  fd(3,1) = x*(y(5)*y(5)-ome2)
  fd(3,2) = x*y(5)*y(5)*rho
  fd(3,3) = -rho - x*rhod
  fd(3,5) = 2.*x*y(5)*y(1) + 2*x*y(5)*rho*y(2)

  f(4) = x*y(5)*y(5)*y(2) + x*y(1)
  fd(4,1) = x
  fd(4,2) = x*y(5)*y(5)
  fd(4,5) = 2.*x*y(5)*y(2)
  
  f(5) = 0.
  
!   if (verbos == 1) then
!      write(*,'(a,2x,i5,11f10.3,e14.4)') '!!!', i, x, ome2, k, rho, rhod, rhodd, psid, c1, a1, a2, c2, a3
!   endif

  
!end_smhr

!      write(*,*) '- Ending rhs - '

end subroutine rhs
      

