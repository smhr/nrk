! *********************************************************************     
subroutine ODEsolve
use initializeNRK
  implicit none

  double precision, dimension(nn) :: x
  double precision, dimension(ii,nn) :: y
  double precision, dimension(ii,3) :: ea
  integer, dimension(ii) :: v
  integer :: i,j,iter, nprevious, n
  double precision :: amean
  
  double precision, dimension(nn) :: rho_val
  double precision :: err
  character(len=40) :: dir_command
  character(len=100) :: sweep_blanks, output_file, not_converged_file, rhs_file
  double precision, dimension(nn) :: bz, dbz
  integer :: verbos
  double precision :: l1, r1, r2, r3
!   double precision, dimension(nn-1) :: dbz
  
  external rhs,bc,am,amd

  verbos = 0
! *** Driver program to perform the solution of ii coupled nonlinear
! *** ODEs in a two-point BVP of the form
! *** Sum_j M_ij d y_j/dx = f_i (x,y)
! *** with ka boundary conditions at x_a and kb boundary conditions at x_b
! *** on any mesh x(1:nn)
  call printPARAMETER ()
  write(*,*) '- Starting ODEsolve - '
  dir_command = 'mkdir brf_results'
  call system ( dir_command )
  dir_command = 'mkdir not_converged'
  call system ( dir_command )
  dir_command = 'mkdir rhs_results'
  call system ( dir_command )
!   call readPARAMETER ()
!   call printPARAMETER ()
  open(21,file='om-k.dat', access='append')

! *** Call to routines to initialize the mesh, and then to input the guess
! print*,'kkkllll'
! open(91,file='rho.dat')
! call mesh(x,rho_val)
!   do i=1,nn
!        write(91,*) x(i), rho_val(i)
!   enddo
! flush(91)

! call density (rho_val)
! print *,'magnetic_field =', magnetic_field
! wave_n = 0.55
wave_n = wave_n_up
! do while ( wave_n < 0.6 .and. wave_n > 0. )
x = 0.d0
nprevious = 2
call mesh(x,nprevious)
!   do while ( wave_n <= wave_n_up .and. wave_n >= wave_n_low )
   do while ( wave_n >= wave_n_up .and. wave_n <= wave_n_low )
!     print*,'pppppppppp'
    ome2 = ome2_up
    omegg2 = - ome2
    omegg = dsqrt(omegg2)
    print*,"in ODE wave_n, ome2, eta =", wave_n, ome2, eta
   do while ( ome2 <= ome2_up .and. ome2 > ome2_low )
!          print*,'fffffffffffff'
         y = 0.d0
         ea = 0.d0; v = 0
         
!          call mesh(x,nprevious)
         call guess(x,y,nprevious)
         ! *** Call to routine to initialize the permutation vector v.

         call initializev(v)

         ! **** Newton-Raphson-Kantorovich loop 
   
         do iter=1,maxiter           !solution of linear equation
!              write(6,'(a20,f8.3)')'in odesolve:,x(2)',x(2)
             call mynrk(x,y,bc,am,amd,ea,v)
             do j=1,ii
                print*, "j, ea",j, ea(j,1)
             enddo

             if(v(1).eq.0) then 
                    write(6,*)'ERROR in V vector.'
                    exit   ! there is an error!
             endif
!              print*,'############'
!              write(6,'(a18,i4,2f8.3)') 'iter, k, ome2 =', iter, wave_n, ome2

             ! ****** Writing out the convergence properties
!              write(*,*) 'Average error at iteration',iter,' is', amean(ea)
         
!              open(33,file='ea.dat',status='unknown')
!              do i=1,ii
!                   write(33,*) (ea(i,j),j=1,3)
!              enddo
!              close(33)
     
             ! ****** If error is small enough, exit loop
!              if(amean(ea).lt.acy) go to 100
!             if (amean(ea).lt.acy) then
!             if (amean(ea).lt.acy .and. dabs((y(6,nn)+y(6,nn-1)/2.d0)) <= 1.d-4 .and. dabs((y(7,nn)+y(7,nn-1)/2.d0)) <= 1.d-4) then
             if (amean(ea).lt.acy) then
!              if (err.lt.acy) then
            write(6,'(a48,3f15.6,i6,e18.5)') 'yes: 1st_k, final_k, ome2, iter ,amean(ea) =', wave_n, dabs(y(5,nn)), ome2, iter &
                                                                                         & ,amean(ea)
                    write(21,'(4f18.10,i5,e18.5)') wave_n, dabs(y(5,nn)), ome2, eta, iter, amean(ea)
                    write(6,'(a27,f18.10)') 'eta in ODEsolve is', eta
                    write(output_file , '( a13, f10.6, a1, f10.6 )' ) 'brf_results/',dabs(ome2),'_',dabs(y(5,nn))
                    flush(21)
                    nprevious = 1
                    call printresult(x,y)
                    output_file = sweep_blanks(output_file)
                    open(44,file=output_file)
                    write(rhs_file , '( a13, f10.6, a1, f10.6 )' ) 'rhs_results/',dabs(ome2),'_',dabs(y(5,nn))
                    rhs_file = sweep_blanks(rhs_file)
                    open(47,file=rhs_file)
                    do i=1, nn
                       write(44,'(f9.4,7f17.12)') x(i), y(1,i), y(2,i), y(3,i), y(4,i), dabs(y(5,i)), y(6,i), y(7,i)
!                       call print_rhs (x(i), y(:,i), i, l1, r1, r2, r3)
!                       write(47,'(f9.4,4f17.12)') x(i), l1, r1, r2, r3
                    enddo
                    close(44); close(47)
!                     stop
                    exit
             endif 
!              nprevious = 1
!              call mesh(x,nprevious)
!              call guess(x,y,nprevious)
           write(6,'(a,2x,3f15.6,i6,e18.5)')'NO: 1st_k, final_k, ome2, iter, amean(ea) =', wave_n, dabs(y(5,nn)), ome2, iter &
                                                                                       & , amean(ea)
             flush(6)
            if (verbos == 1) then
             write(not_converged_file , '( a15, f10.6, a1, f10.6, a1, i3 )' ) 'not_converged/',wave_n,'_',ome2,'_',iter
                    not_converged_file = sweep_blanks(not_converged_file)
                    open(45,file=not_converged_file)
                    do i=1, nn
                       write(45,'(8f17.8)') x(i), y(1,i), y(2,i), y(3,i), y(4,i), dabs(y(5,i)), y(6,i), y(7,i)
                    enddo
            endif
         enddo
        ome2 = ome2 - ome2_step
   enddo
   wave_n = wave_n - wave_n_step
enddo

!!! call printresult(x,y)
close(21)
! *** If it didn't converge no not print results out
write(*,*) '- End ODEsolve - '

end subroutine ODEsolve


! *********************************************************************


 










