subroutine nuclear_wavefkt(adb,chi0)

 use data_grid
 use pot_param
 use data_au

implicit none
  include "/usr/include/fftw3.f"

  integer:: I, J, K, M, V, G   
  integer:: istep  
  integer*8 planF, planB
     
  double precision:: dt2
  double precision:: E, E1, norm
  double precision:: CONS, thresh
  double precision:: dummy
  double precision:: adb(NR, Nstates),chi0(nr,vstates)
 
  double precision, allocatable, dimension(:):: vprop
  double precision, allocatable, dimension(:):: psi, psi1
  double precision, allocatable, dimension(:,:):: ref


 open(99,file='chi0.out',status='unknown')
 open(98,file='Evib.out',status='unknown')

  allocate(psi(NR), psi1(NR))
  allocate(vprop(NR), ref(NR, Vstates))


 call dfftw_plan_r2r_1d(planF, NR, psi, psi, FFTW_R2HC, FFTW_ESTIMATE)
 call dfftw_plan_r2r_1d(planB, NR, psi, psi, FFTW_HC2R, FFTW_ESTIMATE)


  dt2 = dt*5
  thresh = 1.d-16 
  istep = 1e6
 


!............... Main Propagation ........................                   

  print*
  print*, 'Start of energy calculation...'
  print*



Nloop: do J = 1,1! Nstates ! varying the different adiabatic states


    do i = 1, NR ! new vprop
      vprop(i) = dexp(-0.5d0 * dt2 * adb(i,J))
    end do


Vloop:   do V = 1, Vstates ! loop over the vibrational state


      do i = 1, NR  ! symmetry for the startup function
        psi(i) = exp(kappa * (R(I) -RI )**2) +&
     &   (-1.d0)**(V - 1) * exp(-0.5d0 * (R(I) + RI)**2)
      end do
      

      call integ_r(psi, psi, norm)
      psi = psi / sqrt(norm)
      
      E = 0.d0


!.......... Imaginary Time Propagation ........................



    do K = 1, istep

      psi1 = psi        ! storing wave function of iteration step (N - 1)
      E1 = E            ! storing eigenvalue of iteration step (N - 1)


      if (V.gt.1) then   !projecting out the vibrational ground state...
           do G = 1, (V - 1)

            call integ_r(ref(1:NR,G), psi, norm)

            do i = 1, NR
             psi(i) = psi(i) - norm * ref(i, G)
            end do

           end do
      end if


      psi = psi * vprop
      call dfftw_execute(planF)   
      psi = psi * exp((-dt2 * Pr**2) / (2.d0*m_red))   
      call dfftw_execute(planB)   
      psi = psi / dble(Nr)         
      psi = psi * vprop


      call eigenvalue_r(psi, psi1, E, dt2)
      call integ_r(psi, psi, norm)

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
   do i = 1, nr
    write(102,*) i, psi(i)
  end do
stop
        
10  continue



  print*, 'Surface', J, 'Vibrational state', V, E * au2eV
   write(98,*) v, e*au2ev

    do I = 1, NR
      ref(I,V) = psi(I)             ! storing as reference for the next loop
      write(99,*) R(I)*au2a, ref(I,v)+E*au2eV     
    end do
      write(99,*)
   

 end do Vloop               ! end of vibrational states loop


end do Nloop            ! end of surface loop


  chi0 = ref 

  call dfftw_destroy_plan(planF)
  call dfftw_destroy_plan(planB)

  close(99,status='keep')
  close(98,status='keep')

 deallocate(psi, psi1, vprop, ref)

return
end

!_________________ Subroutines______________________________________


subroutine eigenvalue_R(A, B, E, dt2)      
      
use data_grid
 implicit none  
 double precision:: E, e1, e2, norm
 double precision, intent(in):: dt2, A(nr), B(nr)
 
            
  call integ_r(B, B, norm)  
  e1 = norm
  
  call integ_r(A, A, norm)  
  e2 = norm
  
  
  E = (-0.5d0/dt2) * log(e2/e1)
 

return
end subroutine    
  
! ........................................................
                                                          
subroutine integ_r(A, B, C)

use data_grid  
 implicit none  
 integer I 
 double precision,intent(in):: A(Nr), B(Nr)
 double precision C
  
  C = 0.d0
  
  do I = 1, Nr  
   C = C + A(I) * B(I)   
  end do
  
  C = C * dr
  
return  
end subroutine

