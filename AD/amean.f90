! *********************************************************************
double precision function amean(ea)
  use initializeNRK
  implicit none
  
  double precision, dimension(ii,3) :: ea
  double precision :: sum
  integer :: j
  
  sum=0.d0

  do j=1,ii
!         sum=sum+ea(j,2)
     sum=sum+ea(j,1)
  enddo
  
  amean=sum/dble(ii)
 
end function amean

! ***********************************************************************
