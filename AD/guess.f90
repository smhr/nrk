! **************************************************************
subroutine guess(x,y,nprevious)
use initializeNRK
  implicit none

  double precision, dimension(nn), intent(in) :: x
  double precision, dimension(ii,nn), intent(out) :: y
  integer :: i
  integer, intent(in) :: nprevious
  double precision :: dummyx
  double precision, dimension(nn) :: rho, rhod
  double precision :: pd, b
  
  pd = 1.d0 
  b = magnetic_field
  do i=1,nn
         rho(i) = 1.d0/((1.d0 + x(i)*x(i)/8.d0)**2.d0)
         rhod(i) = -2.d0 * (x(i)/4.d0) / ((1.d0 + x(i)*x(i)/8.d0)**3.d0)
  enddo
!   double precision :: k
      
! **************************************************************
! *** Routine in which the guess is entered.
! *** Input parameters:
!     x : the mesh x(1:nn)
!     y : on exit, the guess for every y(i,n)
!
!     The user is prompted if (s)he wants to use a previously
!     computed solution as the next guess. If not, then the user
!     can input the guess by hand in this routine.
! **************************************************************
!   k = 0.5
!   write(*,*) '- Starting guess - ' 
  
!   do i = 1,100
! 	
! 	 print*,i*1., k*x(i), besi0(k*x(i)), 1./(k*k) * (-1.+besi0(k*x(i)))
! 	 
!   enddo
!   stop 
!   write(*,*) 'Use previous result or no? (yes=0,no=1)'
!   read(*,*) nprevious
!   nprevious = 0
print*,"nprevious =", nprevious
  if(nprevious.eq.2) then
  
     open(12,file='result.dat.org',status='old')
     
    do i=1,nn
!         do j=1,ii-1
           read(12,*) dummyx, y(1,i), y(2,i), y(3,i), y(4,i), y(5,i), y(6,i), y(7,i)
           y(5,i) = wave_n
!            y(7,i) = 0.
!            y(3,i) = 0.
!         enddo
!            y(5,i) = wave_n
!            y(2,i) = -y(1,i)/(wave_n*wave_n) ! phi
     enddo
!      y(7,1) = y(7,2)
     close(12)
  elseif(nprevious.eq.1) then
  
     open(13,file='result.dat',status='old')
     
    do i=1,nn
!         do j=1,ii-1
           read(13,*) dummyx, y(1,i), y(2,i), y(3,i), y(4,i), y(5,i), y(6,i), y(7,i)
           y(5,i) = wave_n 
!         enddo
!            y(5,i) = wave_n
!            y(2,i) = -y(1,i)/(wave_n*wave_n) ! phi
     enddo
!      y(7,1) = y(7,2)
     close(13)
     
  else
     
! ****** To be completed by the user : the guess must be entered 
!        for every point, for every function.

!     do j=1,ii
!        do i=1,nn
!           y(j,i) = 1.d0
!        enddo
!     enddo
 
     do i = 1,nn
	y(1,i) = 1.d0 ! f
! 	y(1,i) = 1./(1.+x(i)) ! f
        y(2,i) = -y(1,i)/(wave_n*wave_n) ! phi
        y(3,i) = 0.d0 ! w
        y(4,i) = 0.d0 ! dphi/dx
        y(5,i) = wave_n  ! k
        y(6,i) = (1.d0 - y(5,i)**2.d0/ome2*pd)*b*y(1,i)/rho(i) &
               & - y(5,i)**2./ome2*b*y(2,i)
        y(7,i) = - (1.d0 - y(5,i)**2.d0/ome2*pd)*b*y(1,i)/rho(i)/rho(i) &
               & * rhod(i)
!         print*,i,y(6,i),y(7,i)
     enddo
    
!     write(6,*) 'selected guess values:' 
!     do i = 1, ii
!       write(6,1) i, y(i,1)
!     1 format ('i,y(i,1)', i3, f15.10)
!     enddo
    
  endif
  

!   write(*,*) '- Ending guess - ' 
      

end subroutine guess
