subroutine mynrk(x,y,bc,am,amd,ea,v)
! USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
use initializeNRK
implicit none

 double precision, dimension(nn) ::  x
 double precision, dimension(ii,nn) ::  y
 double precision, dimension(ii,3) :: ea
 integer, dimension(ii) :: v
 double precision :: am,amd


 double precision, dimension(ii,2) :: f
 double precision, dimension(ii,ii,2) :: fd
 double precision, dimension(ii) :: g
 double precision, dimension(ii,ii) :: gs, p, z
 double precision, dimension (ii,ii/2,0:nn-1) :: s
 double precision, dimension(ii,0:nn) :: b
!  double precision, dimension(ii,1) :: sp
 double precision, dimension(ii,ii/2) :: sp

 integer :: i,j,n
 integer :: in,in1,np1,junk
 integer :: ising, routine_index
 integer :: l,ls,lz
 double precision :: sum1,sum2,dmax,test!,aa
 integer :: nmax
 
 double precision, dimension(ii) :: sy

 external rhs,bc,am,amd
 
! subroutine based on the NRK algorithm described by DOG, provides the 
! updated y vectors after 1 iteration of the algorithm, in y, and 
! the corrections average to that vector done in the iteration in ea. 

! ii = total number of equations
! nn = total number of mesh points
! ka = number of bc at first point (entered in module)
! kb = number of bc at second point ka+kb=ii (entered in module)
! ucy = convergence speed (entered in module)
! acy = accuracy required (entered in module)

! x(nn) is the mesh
! y(ii,nn) is the guess for the solutions. 
!          On exit, it has the updated solution
! rhs is the right hand side routine, user-provided
! bc is the boundary conditions routine, user-provided
! am is the derivative of the rhs routine, user-provided
! amd is the derivative of the derivative of rhs routine, user-provided
! ea(ii,3) is the error (contains an average of the error for each 
! variable 
! v is the swapping vector (keeps track of interchanges in the variables 

! z is the current matrix to be zeroed in the sweep algorithm of dimensions 
! of maximum z(ii,ii) (ka is always smaller or equal to ii)
! p is the current matrix to be pivoted in the sweep algorithm, of dimensions
! of p(ii,ii)
! s is the combination of all stored matrices, dimensions are
! s(ii,kb,0:nn+1)
! b is the combination of the rhs vectors b(ii,0:nn+1)

!      write(*,*) ' - Starting mynrk - '

! test solvability
      if(ka+kb.ne.ii) goto 1000
      if(kb.gt.ka) goto 1004
      
     
! empty boundary conditions and derivatives matrices

      do i=1,ii
         do j=1,ii
            gs(i,j) = 0.d0
            fd(i,j,1) = 0.d0
         enddo
      enddo

! set storage indices
      in=1
      in1=2

! initialize error flags
      
      ising=0

! fill up v array 

      do i=1,ii
         if(v(i).eq.0) v(i) = i
!          write(*,*) 'in mynrk: v(',i,')=',v(i)
      enddo
         
! First bc: read in matrices, pivot P_0 and store S_0

      call bc(x(1),x(nn),y(1,1),y(1,nn),g,gs)
      call readinbc1(p,s(1,1,0),b(1,0),g,gs,v)
      routine_index = 0
      call pivot(ka,kb,p,s(1,1,0),b(1,0),ising,routine_index)
      if(ising.eq.1) goto 1001

! Next n blocks: read in blocks, zero Z_n using S_(n-1), update
! P_n and S_n, then pivot P_n storing S_n
!       write(6,'(a10,f8.3)')'1:: x(1)',x(1)
      call rhs(x,y,f,fd,1)
!       write(6,'(a10,2f8.3)')'x(1),x(2)',x(1),x(2)
      do n=1,nn-1

!         write(*,*) ' Current meshpoint ' , n
         np1=n+1
         do l=1,ii
            do j=1,ii
               fd(l,j,in1) = 0.d0
            enddo
         enddo
!           write(*,*) ' Current meshpoint ' , n
! 	 write(6,'(a15,i4,f8.3)')'2:: np1,x(np1)',np1,x(np1)
! 	 write(*,*) '1np1' , np1
         call rhs(x(np1),y(1,np1),f(1,in1),fd(1,1,in1),np1)
!          write(*,*) '2np1' , np1
!          call flush(6)
!          write(6,'(a11,i4,f8.3)')'3:: np1,x(np1)',np1,x(np1)
         call flush(6)
         call readinm(x(n),y(1,n),z,p,s(1,1,n),b(1,n),n,f,fd,in,in1,v,am,amd)

         if(n.eq.1) then
            ls=ka
         else
            ls=ii
         endif
         lz=ii
         call zero(lz,ka,z,p,kb,ls,s(1,1,n-1),b(1,n-1),b(1,n),ii)
         routine_index = 1
         call pivot(ii,kb,p,s(1,1,n),b(1,n),ising,routine_index)
         if(ising.eq.1) goto 1002

         junk=in
         in=in1
         in1=junk
      enddo

! Last block: read Z_(nn), zero it updating P_(nn) then 
! pivot P_(nn) updating b_(nn)

      call readinbc2(z,p,b(1,nn),g,gs,v)
      ls=ii
      lz=kb
      call zero(lz,ka,z,p,kb,ls,s(1,1,nn-1),b(1,nn-1),b(1,nn),ii)
      routine_index = 2
      call pivot(kb,0,p,sp,b(1,nn),ising,routine_index)
      if(ising.eq.1) goto 1003

! Backsusbtitute to obtain the solutions

      call backsub(b,s)

! evaluate error
!  aa = 1.
!  if (wave_n < 0.3) aa = -1.
!  print*,'aa=',aa
      do i=1,ii
         
         do n=1,nn
!             do while (dabs(b(i,n)) > dabs(y(v(i),n)))
! !                y(v(i),n)=y(v(i),n)+ucy/2.*b(i,n)
!                b(i,n) = b(i,n)/2.
!             enddo
!                print*, "i,n,b(i,n)",i, n, b(i,n)
               if (b(i,n) /= b(i,n)) then 
                  print*, "i,n,b(i,n)",i, n, b(i,n)
                  stop
               endif
               y(v(i),n)=y(v(i),n)+ucy*b(i,n)
         enddo
      enddo


      do i=1,ii
         ea(i,1) = 0.d0
!          ea(v(i),1) = 0.d0
         sum1=0.d0
         sum2=0.d0
         dmax=0.d0
         nmax=0
         do n=1,nn
            sum1=sum1+dabs(b(i,n))
            sum2=sum2+dabs(y(v(i),n))
            if(dabs(y(v(i),n)).gt.1.d-15) then
!                test=dabs(b(i,n))/(dabs(y(v(i),n)+1.d-14))
               test=dabs(b(i,n))/(dabs(y(v(i),n)+2.2204460492503131D-016))
            else 
               test = 0.d0
            endif
            if(test.gt.dmax) then
               dmax = test
               nmax = n
            endif
         enddo
         sy(v(i)) = sum2/nn
!          ea(v(i),1) = sum1/sum2
         ea(v(i),1) = sum1/(sum2+2.2204460492503131D-015)
         ea(v(i),2) = dmax
         ea(v(i),3) = nmax 
      enddo

! substitute y for the updated solutions
!       do i=1,ii
!          do n=1,nn
!             y(v(i),n)=y(v(i),n)-ucy*b(i,n)
!          enddo
!       enddo
!       do i=1,ii
!       write (*,*) 'v(i),sy(v(i)),ea(v(i),1)=', v(i), sy(v(i)), ea(v(i),1)
!          do n=1,nn
!             y(v(i),n)=y(v(i),n)+(ucy/max(ucy,ea(v(i),1)))*b(i,n)
!          enddo
!       enddo


!      write(*,*) ' - Ending mynrk -'

      goto 2000

 1000 write(*,*) ' ka + kb is different from I'
      v(1)=0
      goto 2000
 1001 write(*,*) 'First boundary condition matrix singular'
      v(1)=0
      goto 2000
 1002 write(*,*) 'Pivot matrix singular at meshpoint n = ',n
      v(1)=0
      goto 2000
 1003 write(*,*) 'Second boundary condition matrix singular'
      v(1)=0
      goto 2000
 1004 write(*,*) 'kb > ka : reverse mesh and try again'
      v(1)=0
      goto 2000

2000 continue
      
    end subroutine mynrk


! **********************************************************************

subroutine readinbc1(p,s,b,g,gs,v)
use initializeNRK
  implicit none

! subroutine reads in the matrices to be zeroed, pivoted, stored, and the 
! rhs matrices from the user supplied routines rhs, bc, am and amd, 
! in the main mesh sweeps. 
! 
! p = matrix to be pivoted (actual dimensions ka x ka)  
! b = rhs matrix 
! s = matrix to be stored 
! g = bc at meshpoint x(1)
! gs = jacobian of bc at meshpoint x(1)
! v = solutions swapping matrix dimensions are v(ii=ka+kb)

  double precision, dimension(ii,ii) :: p,s,gs
  double precision, dimension(ii) :: b,g
  integer, dimension(ii) :: v
  integer :: j

!      write(*,*) ' - Starting readinbc1 - '
  do j=1,ii
!          write(*,*) 'in readinbc1: v(',j,')=',v(j)
  enddo
      
  do j=1,ka
!      print*,j,gs(1:ka,j),'###',gs(1:ka,v(j)),'***'
     p(1:ka,j) = gs(1:ka,v(j)) 
  enddo
!   print*,gs(1:ka,2)
!   print*,gs(1:ka,5)
  do j=ka+1,ka+kb
     s(1:ka,j-ka)=gs(1:ka,v(j))
  enddo
  b(1:ka) = - g(1:ka)
      
!      write(*,*) ' - Ending readinbc1 - '

end subroutine readinbc1

! ***********************************************************************

subroutine readinbc2(z,p,b,g,gs,v)
use initializeNRK
  implicit none

! subroutine reads in the matrices to be zeroed, pivoted, stored, and the 
! rhs matrices from the user supplied routines rhs, bc, am and amd, 
! in the main mesh sweeps. 
 
! p = matrix to be pivoted of actual dimension (kb,kb)
! z = matrix to be zeroed 
! b = rhs matrix 
! g = bc at meshpoint x(N)
! gs = bc jacobian at meshpoint x(N)
! v = solutions swapping matrix dimensions are v(ii)

  double precision, dimension(ii,ii) :: z,p,gs 
  double precision, dimension(ii) :: b,g
  integer, dimension(ii) :: v
  integer :: i,j

!      write(*,*) ' - Starting readinbc2 - '

  do i=1,kb
     do j=1,ka
        z(i,j) = gs(i+ka,v(j))
     enddo
     do j=ka+1,ka+kb
        p(i,j-ka)=gs(i+ka,v(j))
     enddo
     b(i) = - g(i+ka)
  enddo
  
      
!      write(*,*) ' - Ending readinbc2 - '

end subroutine readinbc2

! ***********************************************************************
subroutine readinm(x,y,z,p,s,b,n,f,fd,in,in1,v,am,amd)
  use initializeNRK
  implicit none

! subroutine reads in the matrices to be zeroed, pivoted, stored, and the 
! rhs matrices from the user supplied routines rhs, bc, am and amd, 
! in the main mesh sweeps. 
! 
! x = mesh points at n and n+1; dimensions are x(2) 
! y = solutions guess matrix at n and n+1; dimensions are y(ii,2)
! z = matrix to be zeroed of dimension (ii,ka)
! p = matrix to be pivoted of dimension (ii,ii)
! s = matrix to be stored of dimension (ii,kb)
! b = rhs matric of dimension (l)
! n = index of the matrices also equivalent to the number of the current 
! mesh point; n varies between 1 and N
! f = rhs at meshpoints n and n+1, dimensions of f are f(ii,2)
! fd = rhs derivatices at meshpoints n and n+1, dimensions 
! of fd are fd(ii,ii,2)
! in = index of meshpoint n (in = 1 or 2)
! in1 = index of meshpoint n+1 (in1 = 1 or 2)
! v = solutions swapping matrix dimensions are v(ii)
! dy = working matrix containing y(i,n+1) - y(i,n), of dimension dy(ii)

  double precision, dimension(2) :: x
  double precision, dimension(ii,2) :: y, f
  double precision, dimension(ii,ii) :: z,p,s
  double precision, dimension(ii) :: b, dy
  double precision, dimension(ii,ii,2) :: fd
  integer, dimension(ii) :: v 
  integer :: n,in,in1
  double precision, external :: am,amd

  double precision :: dx,DYT1,DYT2,DYM
  integer :: i,j,l,np1

! write(*,*) ' - Start readinm - '

! check:

  if(ka+kb.ne.ii) then
     write(*,*) ' problem in readinm' 
     stop
  endif
  
! calculate dy(l)=y(l,n+1)-y(l,n) and dx = x(n+1)-x(n)
  
  np1=n+1

  do l=1,ii
     dy(l) = y(v(l),2)-y(v(l),1)
  enddo

  dx = x(2)-x(1) 

     
! With A_{i,j} given by equations 1-5, we have      
! Z_{i,j} = A_{v(i),v(j)} i=1..ka, j=1..ii
! P_{i,j-ka} = A_{v(i),v(j)}
! P_{i,j+kb} = B_{v(i),v(j)}
! S_{i,j-ka} = B_{v(i),v(j)}
! b_{i} = D_{v(i),v(j)}

  do i=1,ii
     do j=1,ka
        dyT1=0.d0
        dyT2=0.d0
! If you have a very large problem where the LHS is linear, you can comment
! the next loop out to save time.
        do l=1,ii
           dyT1=dyT1+dy(l)*Amd(v(i),v(l),v(j),x(1),y(1,1),n)
           dyT2=dyT2+dy(l)*Amd(v(i),v(l),v(j),x(2),y(1,2),np1)
        enddo
! End comment out. 
        z(i,j)= dx*fd(v(i),v(j),in) + (Am(v(i),v(j),x(1),y(1,1),n) + Am(v(i),v(j),x(2),y(1,2),np1) - dyT1)
        p(i,j+kb) = dx*fd(v(i),v(j),in1)- (Am(v(i),v(j),x(1),y(1,1),n) + Am(v(i),v(j),x(2),y(1,2),np1) + dyT2)
     enddo
     do j= ka+1,ii
        dyT1 = 0.d0
        dyT2=0.d0
! If you have a very large problem where the LHS is linear, you can comment
! the next loop out to save time.
        do l=1,ii
           dyT1=dyT1+dy(l)*Amd(v(i),v(l),v(j),x(1),y(1,1),n)
           dyT2=dyT2+dy(l)*Amd(v(i),v(l),v(j),x(2),y(1,2),np1)
        enddo
! End comment out. 
        p(i,j-ka) = dx*fd(v(i),v(j),in) + (Am(v(i),v(j),x(1),y(1,1),n) + Am(v(i),v(j),x(2),y(1,2),np1) - dyT1)
        s(i,j-ka) = dx*fd(v(i),v(j),in1) - (Am(v(i),v(j),x(1),y(1,1),n) + Am(v(i),v(j),x(2),y(1,2),np1) + dyT2)
     enddo
         
     dyM=0.d0
     do j=1,ii
        dyM=dyM+dy(j)*(Am(v(i),v(j),x(1),y(1,1),n) + Am(v(i),v(j),x(2),y(1,2),np1))
     enddo
     b(i) = - dx*(f(v(i),in)+f(v(i),in1)) + dyM
  enddo
  
  
!      write(*,*) ' - End readinm - '
  
      
end subroutine readinm

! **************************************************************************

subroutine pivot(lp,ks,p,s,b,ising,routine_index)
  use initializeNRK
  implicit none

! Subroutine diagonalises matrix p by implicit partial pivoting, 
! and carries out the same operations on matrix s (to be stored) and
! on the RHS vector b

! Required parameters are the following
! lp = number of lines in matrix to be pivoted
! p = matrix to be pivoted, dimensions are p(lp,lp)
! ks = number of columns in matrix to be stored 
! s = matrix to be stored of dimensions s(lp,kb)
! b = vector to be stored of dimension b(lp)

! Routine written on 02/03/01 by P. Garaud and tested on a simple
! 4x4 matrix in test.f.
  
  double precision, dimension(ii,ii) :: p
  double precision, dimension(ii,ii/2) :: s
  double precision, dimension(ii) :: b, amax
  integer :: lp,ks,ising, routine_index

  integer :: i,j,ip,jp
  double precision :: piv, dum

!      write(*,*) ' - Starting pivot - '
      
! First step, find largest element/row in p
      
  do i=1,lp
     amax(i) = 0.d0
     do j=1,lp
        if(dabs(p(i,j)).gt.amax(i)) amax(i)=dabs(p(i,j))
     enddo
  enddo
  

! Second step, loop over all columns jp the pivoting algorithm

  do jp=1,lp
    
!        write(*,*) 'Pivot number',jp

! 1. Find pivot in column jp (jp = index of column to be pivoted)

     piv=0.d0
     
     do i=jp,lp
        dum= p(i,jp)/amax(i)
        
        if(abs(dum).gt.piv) then
           piv=abs(dum)
           ip = i           ! ip = index of row which holds pivot
        endif
        
     enddo
     
     if(piv.eq.0.d0) goto 1000

!         write(*,*) 'pivot = ', piv, 'in row = ',ip


! 2. Swap row ip with row jp if ip.ne.jp, swap amax, and swap rows in 
! matrices to be stored s and b

     if(ip.ne.jp) then
!            write(*,*) 'swap lines'
       
        do j=jp,lp
           dum = p(ip,j)
           p(ip,j) = p(jp,j)
           p(jp,j) = dum
        enddo
        
        dum = amax(ip)
        amax(ip)=amax(jp)
        amax(jp) = dum
        
        do j=1,ks
           dum = s(ip,j)
           s(ip,j) = s(jp,j)
           s(jp,j) = dum
        enddo
        
        dum = b(ip)
        b(ip)=b(jp)
        b(jp) = dum
     endif
               
! pivot is now placed at position (ip,jp)

     ip=jp

! 3. Normalise row with pivot

     piv=p(ip,jp)
     if(piv.ne.1.d0) then

        do j=jp,lp
           p(ip,j) = p(ip,j)/piv
        enddo
        
        do j=1,ks
           s(ip,j) = s(ip,j)/piv
        enddo
        
        b(ip)=b(ip)/piv
     endif
     
! 4. Eliminate all elements in the column but the pivot, 
! and keep track of the operations on the other columns and stored vectors
! L_i -> L_i - p_(i,jp) L_ip 

     do i= 1,lp
        if(i.ne.jp) then
           if(p(i,jp).ne.0.d0) then
              piv = p(i,jp)
              do j = jp,lp
                 p(i,j) = p(i,j) - piv*p(ip,j)
              enddo
              do j=1,ks
                 s(i,j) = s(i,j) - piv*s(ip,j)
              enddo
              b(i) = b(i) - piv*b(ip)
           endif
        endif
     enddo

! 5. End loop when all columns have been treated that way
 
  enddo

!      write(*,*) ' - Ending pivot - '
  goto 2000

 1000 write(*,*) 'Pivot is null, matrix is singular'
  ising=1
  write(*,*) piv
  do i=1,lp
     do j=1,lp
        write(*,*) 'p(',i,',',j,')=',p(i,j)
     enddo
  enddo
    
2000 continue

! if (routine_index == 2) then
!   do i=1,lp
!      do j=1,lp
!         write(*,*) 'p(',i,',',j,')=',p(i,j)
!      enddo
!   enddo
! endif

end subroutine pivot

! ************************************************************************
subroutine zero(lz,ka,z,p,kb,ls,s,bs,bn,ii)
  implicit none

! Subroutine empties matrix z by substracting with previous iteration's
! result, and carries out the same operations on the first kb 
! columns of matrix p and on the RHS vector b

! Required parameters are the following
! lz = number of lines in matrix to be emptied
! ka = number of columns in matrix to be emptied
! z = matrix to be emptied 
! p = matrix to be pivoted later
! ls = number of lines in stored matrix of previous iterations
! kb = number of columns in stored matrix of previous iterations
! s = stored matrix from previous iteration 
! kb = number of columns in matrix to be stored 
! bs = previous (stored) vector b of dimension b(ls)
! bn = new vector b to be modified of dimension b(lz)
! ii = first dimension of z,p,and s in calling program

! subroutine written on 02/03/01 by P. Garaud
! tested on a simple case with test2.f and test3.f

  double precision, dimension(ii,ii) :: z,p,s
  double precision, dimension(ii) :: bs,bn
  integer :: ka,kb,ls,lz,ii

  integer :: i,jp,j

!      write(*,*) ' - Starting zero - '

  do i=1,lz                    !loop over lines in z & p
     do jp=1,kb                !loop over columns in p
        do j=1,ka              !loop to get all terms cancelled
           p(i,jp) = p(i,jp) - z(i,j)*s(ls-ka+j,jp)
        enddo
     enddo
     do j=1,ka              !loop to get all terms cancelled
        bn(i) = bn(i) - z(i,j)*bs(ls-ka+j)
     enddo
  enddo
!      write(*,*) ' - Ending zero - '

end subroutine zero

! ********************************************************************
subroutine backsub(b,s)
  use initializeNRK
  implicit none

! routine backsubstitutes solution after partial pivoting in nrk algorithm
! and stores the result into b

! b = rhs vector of dimensions b(ii,0:nn)
! s = stored matrices of dimension s(ii,ii,0:nn)

  double precision, dimension(ii,0:nn) :: b
  double precision, dimension(ii,ii/2,0:nn-1) :: s  

  integer :: i,j,n
  double precision :: sum

! if kb lte ka then one can write directly res=b check with Douglas for 
! mesh swapping

! first loop, over last kb variables, in fully diagonalised matrix

!      write(*,*) '- Starting backsub - ' 

  if(kb.gt.ka) goto 1000

  do i=ii,ka+1,-1
     b(i,nn)=b(i-ka,nn)
  enddo

! second loop, sweeping over meshpoints backwards, using the stored matrix
      
  do n=nn-1,1,-1
     do i=ka,1,-1
        sum=0.d0
        do j=1,kb
           sum=sum-s(kb+i,j,n)*b(ka+j,n+1)
        enddo
        b(i,n+1)=b(kb+i,n)+sum
     enddo
     do i=ii,ka+1,-1
        sum=0.d0
        do j=1,kb
           sum=sum-s(i-ka,j,n)*b(ka+j,n+1)
        enddo
        b(i,n) = b(i-ka,n) + sum
     enddo
  enddo
  

! last iterations, for the first boundary conditions
      
  do i=ka,1,-1
     sum=0.d0
     do j=1,kb
        sum=sum-s(i,j,0)*b(ka+j,1)
     enddo
     b(i,1) = b(i,0)+sum
  enddo
  
! b(ii,0) remains unused after backsubstitution

!      write(*,*) ' - Ending backsub - '
  goto 2000

1000  write(*,*) 'kb > ka, cannot use this algorithm to backsubstitute'

2000 continue
end subroutine backsub

! ************************************************************









