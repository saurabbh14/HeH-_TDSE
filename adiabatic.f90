! subroutine to calculate the solution of the electronic Schrödinger equation

subroutine adiabatic_surface(adb, ewf, mu_all)
      
 use data_grid
 use pot_param
 use data_au
 use FFTW3
 use omp_lib

 implicit none
! include "/usr/include/fftw3.f03"
 
  integer I, J, K,L, M, N, G, H, void
  integer*8 planF, planB, istep
  character(100) fn

  real*4 dummy, dummy2, dummy3, dummy4
  double precision dt2
  double precision E, E1, norm, norm1
  double precision CONS,thresh, thresh2
  double precision adb(NR, Nstates)  ! Adiabatic Energy
  double precision ewf(Nx, NR, Nstates), mu(NR,3)
    
  double precision mu_all(NR, Nstates, Nstates)
  double precision, allocatable, dimension(:):: vpropx     
  double precision, allocatable, dimension(:):: psi, psii, psi1  
  double precision, allocatable, dimension(:,:,:):: ref

  open(100,file='HeH+_g0.out',status='replace')
  open(104,file='HeH+_g1.out',status='replace')
  open(102,file='HeH+_bo_potentials.out',status='unknown')
 ! open(103,file='HeH+_g2.out',status='unknown')
  !open(105,file='HeH+_g3.out',status='unknown')
  open(106,file='HeH+_BO.dat',status='unknown')
  open(1000,file='HeH+_trans_dipole.dat',status='unknown')
  open(200,file='HeH+_Ground.dat',status='unknown')
  open(201,file='HeH+_G1.dat',status='unknown')
  !open(202,file='HeH+_G2.dat',status='unknown')
  !open(203,file='HeH+_G3.dat',status='unknown')
  !open(204,file='HeH+_G25.dat')
  
  open(13,file='HeH+_Trans_dipole_all.out',status='unknown')
  open(12,file='dipcurves_3d.dat',status='old')
  read(12,*) dummy, dummy2, dummy3, dummy4
  read(12,*) dummy, dummy2, dummy3, dummy4

  allocate(psi(Nx), psi1(Nx), vpropx(Nx), psii(Nx))
  allocate(ref(Nx, NR, Nstates))

  void =fftw_init_threads( )
  if (void == 0) then
    write(*,*) "Error in fftw_init_threads, quitting"
    stop
  end if

  call fftw_plan_with_nthreads(1)
  call dfftw_plan_r2r_1d(planF, Nx, psi, psi, FFTW_R2HC, FFTW_ESTIMATE)
  call dfftw_plan_r2r_1d(planB, Nx, psi, psi, FFTW_HC2R, FFTW_ESTIMATE)


  dt2 = dt!/10.
  xeq = 2. / au2a

  thresh = 1.d-18 
  thresh2 = 1.d-18 
  istep = 1d6
 
  print*, 'Start of surface scan...'
  print*



states: do N = 1, Nstates       ! scanning the different states
     

   do I = 1, Nx ! choosing the right symmetry for the startup function  
    psi(I) = exp(-10.d0 * (x(I) - xeq)**2) +&
      &  (-1.d0)**(N - 1) * exp(-10.d0* (x(I) + xeq)**2)
   end do
                
   call integ_real(psi, psi, norm)      
      psi = psi / sqrt(norm)

   call integ_real(psi, psi,norm)
   print*, 'norm=', norm   
  print*,'Imaginary Time Propagation, State', N

  !$OMP PARALLEL DO DEFAULT(NONE) FIRSTPRIVATE(vpropx, E, E1,psi1, norm) & 
  !$OMP SHARED(dt2, pot, ewf, thresh, adb,istep, N, planF, planB, px, m_eff, ref) &
  !$OMP FIRSTPRIVATE(psi) 
 
  do I = 1, NR   ! scanning the potential surface -- vary nuclear coordinates   

   do J = 1, Nx 
    vpropx(J) = exp(-0.5d0 * dt2 * pot(I,J))   
   end do
   

 !  psi = psii    ! setting wave function to initial value again
   E = 0.d0     ! setting initial eigenvalue to zero
        

! ........... Imaginary Time Propagation ........................


  do K = 1, istep
  

    psi1 = psi  ! storing wave function of iteration step (N - 1)	    
    E1 = E      ! storing eigenvalue of iteration step (N - 1)


    if (N.gt.1) then     ! projecting out the lower states	   
       do G = 1, (N - 1)

          call integ_real(ref(1:Nx, I, G), psi, norm)
          psi = psi - norm * ref(:, I, G)     
      end do
    end if
 


    psi = psi * vpropx
    call dfftw_execute_r2r(planF,psi,psi)   
    psi = psi * exp((-dt2 * Px**2) / (2.d0*m_eff))   
    call dfftw_execute_r2r(planB,psi,psi)   
    psi = psi / dble(Nx)
    psi = psi * vpropx
  
          
    call eigenvalue_real(psi, psi1, E, dt2) ! calculating the eigenvalue          
    call integ_real(psi, psi, norm) ! normalization of wave function
     
    psi = psi / sqrt(norm)
  
               
    if(abs(E - E1).le.thresh) then     
            goto 10   
    end if

 end do
        
              

   print*,'Iteration not converged!'
   print*,'Program stopped!'
   print*
   print*,'E =', E / cm2au
   print*,'E1 =', E1 / cm2au
   print*,'thresh =', thresh / cm2au
   print*,'step =', K
        
stop
        
10  continue


    adb(I,N) = E
   ! write(N,*) R(I), E
    do J = 1, Nx      
      ref(J,I,N) = psi(J)
    end do
      
                                                                          
   end do                ! end of the R loop
   !$OMP end parallel do 


end do states        ! end of loop over the states




    do I = 1, NR      
      write(106,*) sngl(R(I) *au2a), sngl(adb(I,:)*au2eV)!, sngl(adb(I,2)*au2eV), &
          ! &sngl(adb(i,3)*au2eV), sngl(adb(i,4)*au2eV), ad
    end do


   ewf = ref    ! copy electronic wave- functions into new array
      
            

  mu = 0.d0 ! Dipolemoment

  
      do I = 1, NR    
       do J = 1, Nx
          mu(I,1) = mu(I,1) + ref(J,I,1) * x(J) *ref(J,I,1)
          mu(I,2) = mu(I,2) + ref(J,I,2) * x(J) *ref(J,I,2)
          mu(I,3) = mu(I,3) + ref(J,I,2) * x(J) *ref(J,I,1)! trans. dipole
       end do 
          mu(I,1) = mu(I,1) * dx 
          mu(I,2) = mu(I,2) * dx 
          mu(I,3) = mu(I,3) * dx ! transition dipole
       !read(12,*) dummy, mu(I,1), mu(I,3), mu(I,2)
       write(1000,*) sngl(R(I) *au2a), sngl(mu(I,1)),sngl(mu(I,2)), sngl(mu(i,3))
      end do        
  
      !transition dipole moments of all states
    do L=1,Nstates 
     do M=L+1,Nstates
      write(fn,fmt='(i0,i0,a)') L,M,'.out'
      write(*,*) fn
      open(unit=2000, file=fn, form='formatted')
      do I=1,NR
         mu_all(I,L,M)=sum(ref(:,I,L)*x(:)*ref(:,I,M))*dx
         write(2000,*) R(I)*au2a, mu_all(I,L,M)    
      enddo
      close(2000)
     enddo
    enddo


! wave functions...

   do I = 1, NR  
    !read(106,*) dummy, adb(i,1), adb(i,2), adb(i,3), adb(i,4)
    !read(1000,*) dummy, mu(i,1), mu(i,2), mu(i,3)  
    do J = 1, Nx   
      write(200,*) ref(J,I,1)
      write(201,*) ref(J,I,2)
   !   write(202,*) ref(J,I,3)
     ! write(203,*) ref(J,I,4)
      !write(204,*) ref(J,I,25)
    end do           
   end do   

   ewf = ref
  ! adb = adb/au2eV


   do I = 1, NR, 2 
     write(102,*) sngl(R(I) *au2a), sngl(adb(I,1)*au2eV), sngl(adb(I,2)*au2eV), &
           &sngl(adb(i,3)*au2ev), sngl(adb(i,4)*au2ev)
    do J = nx/4, 3*Nx/4  
      write(100,*) sngl(x(J) *au2a), sngl(R(I) *au2a), sngl(ref(J,I,1))
      write(104,*) sngl(x(J) *au2a), sngl(R(I) *au2a), sngl(ref(J,I,2))
    !  write(103,*) sngl(x(J) *au2a), sngl(R(I) *au2a), sngl(ref(J,I,3))
     ! write(105,*) sngl(x(J) *au2a), sngl(R(I) *au2a), sngl(ref(J,I,4))
    end do       
      write(100,*)
      write(104,*)
     ! write(103,*)
     ! write(105,*)    
   end do   
   print*, R(I)*au2a 

      

  call dfftw_destroy_plan(planF)
  call dfftw_destroy_plan(planB)
 
 deallocate(psi, psii, psi1, vpropx, ref)
      
  close(100,status='keep')
  close(102,status='keep')
  close(103,status='keep')
  close(104,status='keep')
  close(105,status='keep')
  close(106,status='keep')
  close(200,status='keep')
  close(201,status='keep')
  close(202,status='keep')
  close(203,status='keep')
  close(1000,status='keep')
  close(12,status='keep')
return  
end


!_________________ Subroutines______________________________________


subroutine eigenvalue_real(A, B, E, dt2)      
      
use data_grid
 implicit none  
 double precision:: dt2, E, e1, e2
 double precision:: A, B, norm
  
 dimension A(Nx), B(Nx)
            
  call integ_real(B, B, norm)  
  e1 = norm
  
  call integ_real(A, A, norm)  
  e2 = norm
  
  
  E = (-0.5d0/dt2) * log(e2/e1)
  
return
end subroutine    
  
! ........................................................
                                                          
subroutine integ_real(A, B, C)

use data_grid  
 implicit none  
 integer I 
 double precision A(Nx), B(Nx)
 double precision C
  
  C = 0.d0
  
  do I = 1, Nx  
   C = C + A(I) * B(I)   
  end do
  
  C = C * dx
  
return  
end subroutine

