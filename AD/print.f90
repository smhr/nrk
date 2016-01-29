! **********************************************************************
subroutine printresult(x,y)
use initializeNRK
  implicit none

  double precision, dimension(nn) ::  x, rho, rhod
!   double precision, dimension(nn), intent(out) :: bz, dbz
!   double precision, dimension(nn-1), intent(out) :: dbz
  double precision, dimension(ii,nn) :: y
  integer :: i,j

! ***************************************************************
!     Subroutine used only to print the results out. Modify at will
!     but be aware of guess.f where results are read in for next
!     simulation if needed.
! ***************************************************************

!       write(*,*) '- Starting printresult - '

      
  open(50,file='result.dat',status='unknown')
!   open(60,file='f.dat',status='unknown')
!   open(70,file='phi.dat',status='unknown')
!   open(80,file='w.dat',status='unknown')
!   open(90,file='dphi.dat',status='unknown')
!   open(100,file='k.dat',status='unknown')
      
!   do j=1,ii
!      do i=1,nn
!         write(50,*) x(i), y(j,i)
!      enddo
!      write(50,*) ''
!      write(50,*) ''
!   enddo
!   
!   close(50)
  
!   write(*,*)'## 1st k, final k, ome2:', wave_n, y(5,nn), ome2
  do i=1,nn
     write(50,*) x(i), y(1,i), y(2,i), y(3,i), y(4,i), y(5,i), y(6,i), y(7,i)
!      write(60,*) x(i), y(1,i)
!      write(70,*) x(i), y(2,i)
!      write(80,*) x(i), y(3,i)
!      write(90,*) x(i), y(4,i)
!      write(100,*) x(i), y(5,i)
  enddo
!   do i=1,nn-1
! !   do i=1,nn
!        dbz(i) = (bz(i+1)-bz(i))/(x(i+1)-x(i))
! !      write(50,*) x(i), y(1,i), y(2,i), y(3,i), y(4,i), y(5,i), bz(i), dbz(i)
!        write(50,*) x(i), y(1,i), y(2,i), y(3,i), y(4,i), y(5,i), bz(i)
!   enddo
!   write(*,*) '- Ending printresult - '
  close(50)
!   close(60)
!   close(70)
!   close(80)
!   close(90)
!   close(100)
  
end subroutine printresult
         
! *********************************************************************
