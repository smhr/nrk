! **********************************************************************
subroutine print_rhs(x,y,i,l1,r1,r2,r3)
  use initializeNRK
  implicit none
      
  double precision, dimension(ii) :: y, f
  double precision, dimension(ii,ii) :: fd
  double precision :: rho, rhod, rhodd, psid, pd, pdd, b, q1
  double precision :: x, l1, r1, r2, r3, eh
  integer :: i, j, k

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
  
  eh = eta
  if (x > ADlimit) eh = 0.d0
  pd = 1.d0 ; pdd = 0.d0
  b = magnetic_field
!   k = y(5)
  rho = 1.d0/((1.d0 + x*x/8.d0)**2.d0)
  rhod = -2.d0 * (x/4.d0) / ((1.d0 + x*x/8.d0)**3.d0)
  rhodd = 3.d0/8.d0*x*x / ((1.d0 + x*x/8.d0)**4.d0) - 0.5d0/((1.d0 + x*x/8.d0)**3.d0)
  psid = x/(2.d0+x*x/4.d0)
  q1 = -b*b*y(5)*y(5)/(omegg*(eh*b*b*y(5)*y(5)-omegg))

  f = 0.d0; fd = 0.d0
  ! matrix f coefficients
!   c1 = (k*k*b*b)/ome2
!   a1 = x*(-c1/rho*(pdd-2.*pd/rho)*rhod + (pdd-2.*b*b/(rho*rho))*rhod + psid)
!   a2 = x*(c1/rho * rhod)
!   a3 = x*b*b/ome2 * (k*k + rhodd/rho - 2.*(rhod/rho)**2.) - b*b/ome2 * rhod/rho - x*rho
  
  ! matrix fd coefficients
!   c2 = (2.*k*b*b)/ome2
  
!   write(6,'(a25,i8,4f8.3)')'i,x,rho,rhod,psid',i,x,rho,rhod,psid

!!!!!!!!!!!!! start of new f matrix

  f(1) = rho*y(3) - y(1)*pdd*rhod - y(1)*psid + q1*y(3)
  fd(1,1) = - pdd*rhod - psid
  fd(1,3) = rho + q1
  fd(1,5) = (-2.d0*b*b*y(5)*y(3)*(omegg*(eh*b*b*y(5)*y(5)-omegg)) &
          & -2.d0*omegg*eh*b*b*y(5)*(-b*b*y(5)*y(5)*y(3))) &
          & /(omegg*(eh*b*b*y(5)*y(5)-omegg))**2.d0
  
  f(2) = y(4)
  fd(2,4) = 1.d0
   
  f(3) = -y(3)*rho + x*(y(5)*y(5)-ome2)*y(1) + x*y(5)*y(5)*rho*y(2) &
       & - x*y(3)*rhod
  fd(3,1) = x*(y(5)*y(5)-ome2)
  fd(3,2) = x*y(5)*y(5)*rho
  fd(3,3) = -rho - x*rhod
  fd(3,5) = 2.d0*x*y(5)*y(1) + 2*x*y(5)*rho*y(2)

  f(4) = x*y(5)*y(5)*y(2) + x*y(1)
  fd(4,1) = x
  fd(4,2) = x*y(5)*y(5)
  fd(4,5) = 2.d0*x*y(5)*y(2)
  
  f(5) = 0.d0
  
  f(6) = y(7)
  fd(6,7) = 1.d0
  
  f(7) = eh*b*b*y(6)*y(5)*y(5)*x - omegg*y(6)*x &
         & - b*x*(-omegg*y(1)/rho - y(5)*y(5)*y(2)/omegg & 
         & - y(5)*y(5)*y(1)*pd/omegg/rho &
         & + y(3)/omegg/rho*rhod)
  fd(7,1) = -b*x*(-omegg/rho - y(5)*y(5)*pd/omegg/rho)
  fd(7,2) = b*x*y(5)*y(5)/omegg
  fd(7,3) = -b*x*rhod/omegg/rho
  fd(7,5) = 2.d0*eh*b*b*y(6)*x*y(5) - b*x*(-2.d0*y(5)*y(2)/omegg &
            & -2.d0*y(5)*y(1)*pd/omegg/rho)
  fd(7,6) = eh*b*b*y(5)*y(5)*x - omegg*x
  
  !!!!!!!!!!!!!!
  l1 = eh*b*b*y(7)
  r1 = eh*b*b*y(6)*y(5)*y(5)*x
  r2 = omegg*y(6)*x
  r3 = b*x*(-omegg*y(1)/rho - y(5)*y(5)*y(2)/omegg & 
    &              - y(5)*y(5)*y(1)*pd/omegg/rho &
    &              + y(3)/omegg/rho*rhod)
    
    write(*,'(a15,i5,2f17.12)') 'In meshpoint ', i, eh, b
    
write(*,'(11f13.8)') x, ADlimit, b, rho, rhod, rhod/rho, rhodd, psid, q1, y(3), omegg 
!!!!!!!!!!!!! end of new f matrix
!!!!!!!!!!!!! print f and fd matrices
 print*,"++++++++++++++++"
! print*, "omegg2", omegg2
 do j = 1, ii
    write (*,'(a35,i3,f14.4,5e11.3)') &
    &     "x, j, y(j), rho, rhod, psid, f(j)", &
    &      j, x, y(j), rho, rhod, psid, f(j)
 enddo
 print*,"================"
 do j = 1, ii
    write (*,'(7e12.3)') (fd(j,k), k = 1, ii)
 enddo
!!!!!!!!!!!!!
! !   f(1) = y(3)*rho - y(1)*psid
!   f(1) = -a1*y(1) - a2*y(2) - a3*y(3)
! !   fd(1,1) = -psid
!   fd(1,1) = -a1
!   fd(1,2) = -a2
! !   fd(1,3) = rho
!   fd(1,3) = -a3
!   fd(1,5) = x*(c2/rho*(pdd-2.*pd/rho)*rhod * y(1) - (c2/rho * rhod) * y(2) - c2 * y(3))
!   
!   f(2) = y(4)
!   fd(2,4) = 1.
!    
!   f(3) = -y(3)*rho + x*(y(5)*y(5)-ome2)*y(1) + x*y(5)*y(5)*rho*y(2) - x*y(3)*rhod
!   fd(3,1) = x*(y(5)*y(5)-ome2)
!   fd(3,2) = x*y(5)*y(5)*rho
!   fd(3,3) = -rho - x*rhod
!   fd(3,5) = 2.*x*y(5)*y(1) + 2*x*y(5)*rho*y(2)
! 
!   f(4) = x*y(5)*y(5)*y(2) + x*y(1)
!   fd(4,1) = x
!   fd(4,2) = x*y(5)*y(5)
!   fd(4,5) = 2.*x*y(5)*y(2)
!   
!   f(5) = 0.
  
!   if (verbos == 1) then
!      write(*,'(a,2x,i5,11f10.3,e14.4)') '!!!', i, x, ome2, k, rho, rhod, rhodd, psid, c1, a1, a2, c2, a3
!   endif

  
!end_smhr

!      write(*,*) '- Ending rhs - '

end subroutine print_rhs
      

