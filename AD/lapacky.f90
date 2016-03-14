subroutine lapacky (x,y,bc,am,amd,ea,v)

use initializeNRK
! USE LA_PRECISION, ONLY: WP => DP
! USE F95_LAPACK, ONLY: LA_GESV
implicit none

 double precision, dimension(nn) ::  x
 double precision, dimension(ii,nn) ::  y
 double precision, dimension(ii,3) :: ea
 integer, dimension(ii) :: v
 double precision :: am,amd
 
 double precision, dimension(ii,ii) :: AMAT, BMAT, gs
 double precision, dimension(ii) :: DVEC, dy, g
 double precision, dimension(ii,2) :: f
 double precision, dimension(ii,ii,2) :: fd
 double precision, dimension(ii*nn) :: b
!  REAL(WP), allocatable :: b(:)
 double precision, dimension(ii*nn,ii*nn) :: wMAT
!  REAL(WP), dimension(ii*nn,ii*nn) :: wMAT
!  REAL(WP), dimension(ii*nn,1) :: bMKL
 
 integer :: nstore,i,j,l,ll,n
 integer :: in,in1,junk
 integer :: np1
 double precision :: dx,dyT1,dyT2,DYM,am1,am2
 double precision :: sum1,sum2,dmax,test
 integer :: nmax
 integer :: i1, j1
 
 integer :: info, lda, ldb, nrhs
 integer, dimension(ii*nn) :: ipiv
  
 external rhs,bc,am,amd
 
! test solvability
 if (ka+kb.ne.ii) then
    write(*,*) ' ka + kb is different from I'
    stop
 endif
 
 wMAT = 0.d0
 do i = 1, ii*nn
    b(i) = 0.d0
 enddo
! set storage indices
 in=1; in1=2     
! empty boundary conditions and derivatives matrices
 do i=1,ii
    g(i) = 0.d0
    do j=1,ii
       gs(i,j) = 0.d0
       fd(i,j,in) = 0.d0
    enddo
 enddo
  

! initialize error flags 
!       ising = 0
      
! fill up v array 
 do i = 1, ii
    if (v(i).eq.0) v(i) = i
 enddo
      
! Read in first and last boundary condition matrices and compute the 1st and last block
! in W matrix and in b vector.
 call bc(x(1),x(nn),y(:,1),y(:,nn),g,gs)
!  
 do i = 1, ka
    b(i) = - g(i)
    print*,"i, b(i)",i, b(i)
 enddo
 do j = 1, ii
    do i = 1, ka
!        print*, i, j, v(j), gs(i,v(j))
       wMAT(i,j) = gs(i,v(j))  
!        wMAT(i,j) = 1.0
    enddo
 enddo
 
 do j = 1, kb
     i1 = (nn - 1) * ii + ka + j
     b(i1) = - g(j+ka)
     print*,"i1, b(i1)",i1, b(i1)
 enddo
 do j = 1, ii
       j1 = (nn - 1) * ii + j
       do i = 1, kb
          i1 = (nn - 1) * ii + ka + i
          wMAT(i1,j1) = gs(ka + i,v(j))
       enddo
 enddo
  
!!!!!!!!!!!!!!!!
      
 call rhs(x,y,f,fd,1)
! compute matrices A and B      
      
!  wMAT = 0.d0
! read big matrix wMAT
 do 100 n = 1, nn - 1
    print*,"n =", n
 
    np1 = n + 1
    do l=1,ii
       do j=1,ii
          fd(j,l,in1) = 0.d0
       enddo
    enddo
! calculate dy and dx      
    do l=1,ii
       dy(l) = y(v(l), n + 1) - y(v(l), n)
    enddo
    dx = x(n + 1) - x(n)
    call rhs(x(np1),y(1,np1),f(1,in1),fd(1,1,in1),np1)     
    
    do j = 1,ii
        do i = 1,ii
           dyT1=0.d0
           dyT2=0.d0
!            if(linearlhs.eq.0) then
              do l=1,ii
                 dyT1=dyT1+dy(l)*Amd(v(i),v(l),v(j),x(n),y(1,n),n)
                 dyT2=dyT2+dy(l)*Amd(v(i),v(l),v(j),x(np1),y(1,np1),np1)
!                  write(*,'(a30,3x,4i4,f16.8)') "i, l, j, n,  amd", i, l, j, n, amd(v(i),v(l),v(j),x(1),y(1,1),n)
!                  write(*,'(a30,3x,4i4,f16.8)') "i, l, j, np1, amd", i, l, j, np1, amd(v(i),v(l),v(j),x(2),y(1,2),np1)
              enddo
!            endif
           am1 = Am(v(i),v(j),x(n),y(1,n),n) 
           am2 = Am(v(i),v(j),x(np1),y(1,np1),np1) 
           AMAT(i,j)= dx*fd(v(i),v(j),in) + (Am1 + Am2 - dyT1)
           BMAT(i,j) = dx*fd(v(i),v(j),in1) - (Am1 + Am2 + dyT2)
        enddo
     enddo
         
         do j = 1, ii
            j1 = (n - 1) * ii + j
            do i = 1, ii
               i1 = (n - 1) * ii + i + ka
               wMAT(i1,j1) = AMAT(i,j)
            enddo
         enddo
         do j = 1, ii
            j1 = n * ii + j
            do i = 1, ii
               i1 = (n - 1) * ii + i + ka
               wMAT(i1,j1) = BMAT(i,j)
            enddo
         enddo
         
         do i = 1, ii
            dyM = 0.d0
            do j = 1, ii
               dyM=dyM+dy(j)*(Am(v(i),v(j),x(n),y(1,n),n) + Am(v(i),v(j),x(np1),y(1,np1),np1))
!                write(*,'(a30,3x,3i4,f16.8)') "n, i, j, dym", n, i, j, dym
!                write(*,'(a30,3x,3i4,f16.8)') "i, l, j, n, Am", i, j, n, Am(v(i),v(j),x(n),y(1,n),n)
!                write(*,'(a30,3x,3i4,f16.8)') "i, l, j, np1, Am", i, j, np1, Am(v(i),v(j),x(np1),y(1,np1),np1)
            enddo
            DVEC(i) = dx*(f(v(i),in)+f(v(i),in1)) - dyM
!             write(*,'(a35,3x,3i3,4f16.8)') "i, in, in1, dx, fn, fn1, DVEC(i)", i, in, in1, dx, f(v(i),in), f(v(i),in1), DVEC(i)
         enddo
         do i = 1, ii
            i1 = (n - 1) * ii + i + ka
            b(i1) = -DVEC(i)
            write(*,'(a30,3x,2i4,f16.8)') "n, i1, b(i1)", n, i1, b(i1)
         enddo
         junk=in
         in=in1
         in1=junk
  100 enddo    
  
!  do i = 1, ii*nn
!     bMKL(i,1) = b(i)
!  enddo
! call lapack driver routine
!  CALL LA_GESV(  wMAT, bMKL )
 nrhs = 1 ! number of right hand sides in b
 lda = ii*nn  ! leading dimension of a
 ldb = ii*nn  ! leading dimension of b

    call dgesv(ii*nn, nrhs, wMAT, lda, ipiv, b, ldb, info)
    print*,"info =", info
    print*,"b(1), b(8), b(nn) =", b(1), b(8), b(nn)
!     stop
!  do i = 1, ii*nn
!     b(i) = bMKL(i,1)
!  enddo
!!!!!!!!!!!!!!!!!!
! update solution
 do n = 0, nn - 1 
    do i = 1, ii
!        if (b(i,n) /= b(i,n)) then 
!           print*, "i,n,b(i,n)",i, n, b(i,n)
!           stop
!        endif
       y(v(i),n+1) = y(v(i),n+1) + ucy * b(n*ii + i)
     enddo
 enddo
 
! calculate relative error
    do i = 1, ii
         ea(i,1) = 0.d0
!          ea(v(i),1) = 0.d0
         sum1=0.d0
         sum2=0.d0
         dmax=0.d0
         nmax=0
         do n = 0, nn -1
            sum1=sum1+dabs(b(n*ii + i))
            sum2=sum2+dabs(y(v(i),n+1))
            if(dabs(y(v(i),n+1)).gt.1.d-15) then
!                test=dabs(b(i,n))/(dabs(y(v(i),n)+1.d-14))
               test=dabs(b(n*ii + i))/(dabs(y(v(i),n+1)+2.2204460492503131D-016))
            else 
               test = 0.d0
            endif
            if(test.gt.dmax) then
               dmax = test
               nmax = n
            endif
         enddo
!          ea(v(i),1) = sum1/sum2
         ea(v(i),1) = sum1/(sum2+2.2204460492503131D-015)
         ea(v(i),2) = dmax
         ea(v(i),3) = nmax 
      enddo
      print*,"sum1, sum2", sum1, sum2
!       stop

end subroutine lapacky