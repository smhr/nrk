! **************************************************************
subroutine guess(x,y,nprevious)
use initializeNRK
  implicit none

  double precision, dimension(nn) :: x
  double precision, dimension(ii,nn) :: y
  integer :: i
  integer, intent(in) :: nprevious
  double precision :: dummyx
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

  if(nprevious.eq.1) then
  
     open(12,file='result.dat',status='old')
     
    do i=1,nn-1
!         do j=1,ii-1
           read(12,*) dummyx, y(1,i), y(2,i), y(3,i), y(4,i)!, y(5,i)
           y(5,i) = wave_n 
!         enddo
!            y(5,i) = wave_n
!            y(2,i) = -y(1,i)/(wave_n*wave_n) ! phi
     enddo
     close(12)
     
  else
     
! ****** To be completed by the user : the guess must be entered 
!        for every point, for every function.

!     do j=1,ii
!        do i=1,nn
!           y(j,i) = 1.d0
!        enddo
!     enddo
 
     do i = 1,nn
! 	y(1,i) = 1. ! f
	y(1,i) = 1./(1.+x(i)) ! f
        y(2,i) = -y(1,i)/(wave_n*wave_n) ! phi
        y(3,i) = 0. ! w
        y(4,i) = 0. ! dphi/dx
        y(5,i) = wave_n  ! k
     enddo
    
!     write(6,*) 'selected guess values:' 
!     do i = 1, ii
!       write(6,1) i, y(i,1)
!     1 format ('i,y(i,1)', i3, f8.3)
!     enddo
    
  endif
  

!   write(*,*) '- Ending guess - ' 
      

end subroutine guess
