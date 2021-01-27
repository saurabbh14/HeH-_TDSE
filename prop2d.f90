subroutine propagation_2D(ewf,chi0)

use data_grid
use pot_param
use data_au
use FFTW3
use omp_lib

 implicit none
! include "/usr/include/fftw3.f"

 integer I, J, K
 integer L, M, N, void
 integer*8 planF, planB, planFd 
 integer eR

 double precision dt2, time, deR
 double precision E2, E1, E21, E22
 double precision norm
 double precision norm_outx, norm_outR

 double precision evR, evx, epx, epr
 double precision x_out, R_out
 double precision momt_x,momt_R

 double precision kl, lam !Laser coupling (k*x + lambda*R)*Electric field
 
 double precision, intent(in):: ewf(Nx, NR, Nstates),chi0(nr,vstates)
 double precision vib_pop(vstates)
 double precision:: cof(Nx), cofR(NR)
 complex*16:: corrf, normp, normn
 complex*16, allocatable, dimension(:,:):: psi, kprop, psi0
 double precision, allocatable, dimension(:):: idensR, idensx
 double precision, allocatable, dimension(:):: idenspR, idenspx
 complex*16, allocatable, dimension(:,:):: psi_out_x, psi_out_R
 complex*16, allocatable, dimension(:,:):: psi_dum 
 complex*16, allocatable, dimension(:,:):: psi_out_x1, psi_out_R1


 open(98,file='cof_2d.out',status='unknown')
 open(99,file='cofR_2d.out',status='unknown')
 open(100,file='psi0_2d.out',status='unknown')
 open(200,file='dens_R.out',status='replace')
 open(201,file='dens_x.out',status='replace')
 open(202,file='Pdens_R.out',status='replace')
 open(203,file='Pdens_x.out',status='replace')


 open(500,file='pop1_2d.out',status='replace')
 open(501,file='pop2_2d.out',status='replace')
 open(502,file='pl_pop1_2d.out',status='replace')
 open(503,file='neg_pop2_2d.out',status='replace')
 open(504,file='densR1_2d.out',status='replace')
 open(505,file='densR2_2d.out',status='replace')
 open(506,file='vib_pop_2d.out',status='replace')
 open(600,file='ampl1_2d.out',status='replace')
 open(601,file='ampl2_2d.out',status='replace')
 open(800,file='R_2d.out',status='replace')
 open(801,file='x_2d.out',status='replace')
 open(802,file='norm_2d.out',status='replace')
 open(803,file='corr_2d.out',status='replace')

 open(909,file='field_2d.out',status='replace')
 open(907,file='field1_2d.out',status='replace')
 open(908,file='field2_2d.out',status='replace') 
 open(910,file='pulse.out',status='replace')

 open(1000,file='test.out')
 open(1111,file='accumulation_x.out',status='unknown')
 open(1112,file='accumulation_R.out',status='unknown')
 open(1113,file='nclr_momt_2d.out',status='unknown')
 open(1114,file='elec_momt_2d.out',status='unknown')
 open(1115,file='KER_R2d.out',status='unknown')

 allocate(psi(NR, Nx), idensR(NR), idensx(Nx))
 allocate(kprop(NR,Nx), psi0(nr,nx))
 allocate(idenspr(nr), idenspx(nx))
 allocate(psi_out_x(NR,Nx),psi_out_R(NR,Nx))
 allocate(psi_dum(NR,Nx))
 allocate(psi_out_x1(NR,Nx),psi_out_R1(NR,Nx))


 print*
 print*,'Tuning FFTW...'

 void=fftw_init_threads()
 if (void .eq. 0) then
      print*, 'Error in fftw_init_threads, quitting'
 endif

 call fftw_plan_with_nthreads(4)
 call dfftw_plan_dft_2d(planF, NR, Nx, psi, psi, FFTW_FORWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planB, NR, Nx, psi, psi, FFTW_BACKWARD,FFTW_MEASURE)
 call dfftw_plan_dft_2d(planFd, NR, Nx, psi_dum, psi_dum, FFTW_FORWARD,FFTW_MEASURE)

 print*,'Done.'
 print*


 do I = 1, NR! pR**2 /2* red_mass
   do J = 1, Nx ! px**2 / 2* m_eff 
    kprop(I,J) = exp(-im *dt *((PR(I)**2/(2.d0*m_red))+(Px(J)**2/(2.d0 *m_eff))))
   end do
 end do
 
 do I=1,NR
   do J=1,Nx
   psi(I,J) = ewf(J,I,1)*(0.55d0*chi0(I,1)+0.23d0*chi0(I,2)+0.11d0*chi0(I,3) &
             & + 0.07d0 * chi0(I,4) + 0.04d0*chi0(I,5)+ 2.0450413213653231E-002*chi0(I,6) &
             & + 1.4033977605843639E-002 *chi0(I,7)+ 1.0613089628800979E-002*chi0(I,8) &
             & + 8.8496818475779886E-003*chi0(I,9)+ 8.0611229210941892E-003*chi0(I,10))
!  psi(I,J)= ewf(J,I,1)*exp(kappa*(r(I)-ri)**2)!chi0(i,3)+chi0(i,4)+chi0(i,5) 
   enddo
 enddo


 call cofs(cof, cofR)

 call integ_2d(psi, norm)
 print*,'norm1 =', sngl(norm)

 psi = psi / sqrt(norm)

 call integ_2d(psi, norm)
 print*,'norm2 =', sngl(norm)

 do I = 1, NR/2
   do J = 1, Nx, 4
    write(100,*) sngl(R(I) *au2a), sngl(x(J) *au2a),sngl(abs(psi(I,J))**2)!density?
   end do
    write(100,*)
 end do

 psi0 = psi
 
!______________________________________________________________________
!
!                   Propagation Loop
!______________________________________________________________________


 print*
 print*,'2D propagation...'
 print*

 !psi_out initiation___!
 do J=1,nx
    psi_out_x(:,J)=psi(:,J)*(1.0d0-cof(J))
 enddo
 psi_dum=psi_out_x
 call dfftw_execute(planFd)
 psi_out_x=psi_dum
 do I=1,NR
    psi_out_R(I,:)=psi(I,:)*(1.0d0-cofR(I))
 enddo
 psi_dum=psi_out_R
 call dfftw_execute(planFd)
 psi_out_R=psi_dum
 ! Done________________!

 !Laser coupling parameters____________!
E22=0.0d0
kl=(mn+2.0d0)/(mn+1)
lam=(m2-m1)/mn



 !starting timeloop____!
timeloop: do K = 1, Nt

    time = K * dt

    if(mod(K,200).eq.0) then
      print*,'time:', sngl(time *au2fs)
    end if

   evR = 0.d0
   evx = 0.d0
   epR = 0.d0
   epx = 0.d0
   
    E21 = E01 *exp(-fwhm * (time - t_start)**2)* (cos(omega1 * time+phi1))! &
      ! & - ((2.d0 * fwhm) / omega1) * (time - t_start) * sin(omega1 * time+phi1))
    
!    E22 = E02 *exp(-fwhm * (time - t_start)**2)* (cos(omega2 * time+phi2))! &
      ! & - ((2.d0 * fwhm) / omega2) * (time - t_start) * sin(omega2 * time+phi2))

    
    E2=E21+E22  

 !........................
 

   
    do I = 1, Nr
      do J = 1, Nx
        psi(i,j) = psi(i,j) * exp(-0.5d0*im * dt * (pot(i,j)  &
               & + (kl*x(j)+lam*R(I)) *E2))
      end do
    end do

      call dfftw_execute(planF, psi,psi)
 
      psi = psi * kprop
      psi = psi / sqrt(dble(NR * Nx))

      
        call density(psi,idenspR,idenspx)
        do I = 1, NR
         epR = epR + dble(pR(i) *idenspR(I))
        end do
        epR = epR * dR/ norm

        do J = 1, Nx
         epx = epx + dble(px(j) *idenspx(J))
        end do
         epx = epx * dx/ norm
       
      call dfftw_execute(planB, psi, psi)
       psi = psi / sqrt(dble(NR * Nx))
    

    do I = 1, Nr
      do J = 1, Nx
      psi(i,j) = psi(i,j) * exp(-0.5d0* im * dt * (pot(i,j) + &
               & + (kl*x(j)+lam*R(I)) *E2))
      end do
    end do
 
! ............................................



  call density(psi,idensR,idensx)
  call integ_2D(psi, norm)


   do I = 1, NR
    evR = evR + dble(R(I) *idensR(I))
   end do
    evR = evR * dR

   do J = 1, Nx
    evx = evx + dble(x(j) *idensx(J))
   end do
    evx = evx * dx

    evR = evR / norm
    evx = evx / norm


 call overlap_2D(psi0, psi, corrf)
 call pop_analysis(psi, time, ewf, K)
 call localization(psi, time, ewf, K)
 
    write(506,*) sngl(time*au2fs), (vib_pop(1)), (vib_pop(2)), (vib_pop(3)), (vib_pop(4)),&
              &  vib_pop(5), vib_pop(6), vib_pop(7), vib_pop(8), vib_pop(9)!,&
            ! &  psi_chi(10), psi_chi(11), psi_chi(12)   


    write(800,*) sngl(time *au2fs), sngl(evR *au2a), sngl(epR)
    write(801,*) sngl(time *au2fs), sngl(evx *au2a), sngl(epx)
    write(802,*) sngl(time *au2fs), sngl(norm) 
    write(803,*) sngl(time *au2fs), sngl(abs(corrf))
    write(909,*) sngl(time *au2fs), sngl(E2)
    write(907,*) sngl(time *au2fs), sngl(E21)
    write(908,*) sngl(time *au2fs), sngl(E22)
    
 if(mod(K,50).eq.0) then
  do I = 1, NR/2
    write(200,*) sngl(time *au2fs), sngl(R(I) *au2a),sngl(idensR(I))
  end do
   write(200,*)
  
   do I = NR/2+1, NR
    write(202,*) sngl(time *au2fs), sngl(Pr(i)),sngl(idenspR(I))
   end do
   do I = 1, NR/2
    write(202,*) sngl(time *au2fs), sngl(pr(i)),sngl(idenspR(I))
   end do
   write(202,*)
end if
if(mod(K,10).eq.0) then
  do J = 7*nx/16, 9*Nx/16
    write(201,*) sngl(time *au2fs), sngl(x(j) *au2a), sngl(idensx(J))
  end do
    write(201,*)

!  do j = 9*Nx/16, Nx
!    write(203,*) sngl(time *au2fs), sngl(Px(j)),sngl(idenspx(j))
!  end do
!  do j = 1,7*Nx/16
!    write(203,*) sngl(time *au2fs), sngl(px(j)),sngl(idenspx(j))
!  end do
!   write(203,*)
 end if


 


 !...............CONTINUUM....................................
 ! Electronic coordinate wavefunction prpopagation____!

   do I=1,NR
    do J=1,Nx
      psi_out_x(i,j) = psi_out_x(i,j) * exp(-im * dt *  &
                   &  (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x(j) *E2)
    enddo
   enddo
   psi_out_x=psi_out_x*kprop

 ! Nuclear coordinate wavefunction propagation________!
   do I=1,NR
    do J=1,Nx
     psi_out_R(i,j) = psi_out_R(i,j) * exp(-im * dt *  &
                   &  (1.d0 + 1.d0/(2.d0*m_red +1.d0)) *x(j) *E2)
    enddo
   enddo
    psi_out_R=psi_out_R*kprop
 ! applying cutoff____________________________________!
  do J = 1, Nx
   do I = 1, NR
   psi_out_x1(i,j)=psi(i,j)*(1.0d0-cof(j))*(1.0d0-cofR(i))
   psi_out_R1(i,j)=psi(i,j)*(1.0d0-cofR(i))
   psi(i,j) = psi(i,j) * cof(j) * cofR(i)
   end do
  end do

 ! Fourier transform of the cutoffs____________________!
  psi_dum=psi_out_x1

  call dfftw_execute(PlanFd, psi, psi)

  psi_out_x1=psi_dum
  psi_out_x1=psi_out_x1/sqrt(dble(NR*Nx))
  call integ_2D(psi_out_x1, norm_outx)

  psi_dum=psi_out_R1
  call dfftw_execute(PlanFd, psi, psi)
 
  psi_out_R1=psi_dum
  psi_out_R1=psi_out_R1/sqrt(dble(NR*Nx))
  call integ_2D(psi_out_R1, norm_outR)

! Accumulation________________________________________!
  psi_out_x=psi_out_x+psi_out_x1
  call integ_2D(psi_out_x, norm_outx)

  psi_out_R=psi_out_R+psi_out_R1
  call integ_2D(psi_out_R, norm_outR)

 ! Analysis___________________________________________!
 ! Electronic Momentum
  momt_x=0.0d0
  call density(psi_out_x, idensR, idensx)
  do J=1,Nx
    momt_x=momt_x+dble(px(J)*idensx(J))
  enddo
    momt_x=momt_x*dx
    momt_x=momt_x/norm_outx
    write(1114,*) sngl(time*au2fs), sngl(momt_x), sngl(norm_outx)

 ! Nuclear Momentum   
  momt_R=0.0d0
  call density(psi_out_R, idensR, idensx)
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
    momt_R=momt_R/norm_outR
    write(1113,*) sngl(time*au2fs), sngl(momt_R), sngl(norm_outR)


end do timeloop

 !normalization of outgoing wavefunctions
  call integ_2D(psi_out_x,norm_outx)
  psi_out_x=psi_out_x/sqrt(norm_outx)
  call integ_2D(psi_out_R, norm_outR)
  psi_out_R=psi_out_R/sqrt(norm_outR)

 ! Momemtum distribution
  call density(psi_out_x, idensR, idensx)
  do J=Nx/2+1, Nx
   write(1111,*) sngl(px(J)), sngl(dble(idensx(J)))
  enddo
  do J=1,Nx/2
   write(1111,*) sngl(px(J)), sngl(dble(idensx(J)))
  enddo

  call density(psi_out_R, idensR, idensx)
  do I=NR/2+1, NR
   write(1112,*) sngl(pR(I)), sngl(dble(idensR(I)))
  enddo
  do I=1,NR/2
   write(1112,*) sngl(pR(I)), sngl(dble(idensR(I)))
  enddo
   do I=1, NR/4
    write(1115,*) (pr(I)**2)*au2ev/(2.00d0*m_red), dble(idensR(I)), dble(idensR(I))
   enddo


  momt_R=0.0d0 
  do I=1,NR
    momt_R=momt_R+dble(pR(I)*idensR(I))
  enddo
    momt_R=momt_R*dR
!    momt_R=momt_R/norm_outR
  print*, 'average momentum=', sngl(momt_R), 'a.u.'
  print*, 'average velocity=', sngl((momt_R*au2a)/(mass*au2fs)), 'A°/fs'


do I= -Nt, Nt
    time=I *dt
    E21 = E01 *exp(-fwhm * (time - t_start)**2)* (cos(omega1 * time+phi1))! &
       ! & - ((2.d0 * fwhm) / omega1) * (time - t_start) * sin(omega1 * time+phi1))

    E22 = E02 *exp(-fwhm * (time - t_start)**2)* (cos(omega2 * time+phi2))! &
       ! & - ((2.d0 * fwhm) / omega2) * (time - t_start) * sin(omega2 * time+phi2))

    E2=E21+E22
    write(910,*) sngl(time *au2fs), sngl(E2)
enddo

 deallocate(idensR, idensx, psi, kprop, psi0, idenspx, idenspR)
 deallocate(psi_out_R, psi_out_x, psi_out_R1, psi_out_x1) 



! close(98, status='keep')
 close(100, status='keep')
 close(200, status='keep')
 close(201, status='keep')
 close(202, status='keep')
 close(203, status='keep')
 close(500, status='keep')
 close(501, status='keep')
 close(502, status='keep')
 close(503, status='keep')
 close(550, status='keep')
 close(600, status='keep')
 close(601, status='keep')
 close(800, status='keep')
 close(801, status='keep')
 close(802, status='keep')
 close(803, status='keep')
 close(909, status='keep')
 close(907, status='keep') 
 close(908, status='keep')
 close(910, status='keep')


 call dfftw_destroy_plan(planF)
 call dfftw_destroy_plan(planB)

 

return
end subroutine

!_________________________________________________________

subroutine integ_2D(psi, norm)

use data_grid
 implicit none
 integer I, J

 double precision,intent(out):: norm
 complex*16,intent(in):: psi(NR, Nx)

 norm = 0.d0
 
  do I = 1, NR
   do J = 1, Nx    
      norm = norm + abs(psi(i,j))**2
   end do
  end do

   norm = norm* dx*dR
return
end subroutine


! ..................................................................

subroutine density(psi,idensR,idensx)

use data_grid
use pot_param,only:R0,x0

 implicit none
 integer:: I, J
 double precision,intent(out)::idensx(Nx), idensR(NR) 
 complex*16,intent(in):: psi(NR,Nx)
 
   idensR = 0.d0
   idensx = 0.d0

  do I = 1, NR
    do J = 1, Nx
      idensR(I) = idensR(I) + abs(psi(I,J))**2
    end do
      idensR(I) = idensR(I) *dx
  end do


  do J = 1, Nx
    do I = 1, NR
      idensx(J) = idensx(J) + abs(psi(I,J))**2
    end do
      idensx(J) = idensx(J) *dR
   end do

return
end subroutine

!...................................................

subroutine overlap_2d(psi1, psi2, C)

use data_grid
 implicit none
 integer:: I, J
 complex*16:: C, F(NR)
 complex*16,intent(in):: psi1(NR,Nx), psi2(NR,Nx)

 C = (0.d0, 0.d0); F = (0.d0, 0.d0)

 do I = 1, NR
  do J = 1, Nx

   f(i) = f(i) + conjg(psi1(i,j)) * psi2(i,j)

  end do
 end do

  f = f * dx

  do i = 1, NR
   C = C + F(I)
  end do

  C = C * dR


return
end subroutine

!........................................................................

subroutine pop_analysis(psi, time, ewf, K)

use data_grid
use data_au
use pot_param

 implicit none

 integer:: I, J, K, N

 double precision,intent(in):: time, ewf(Nx, NR, Nstates)
 double precision:: B(Nstates), ax(2), cx(nr,2)
 complex(kind=kind(0.d0)),intent(in):: psi(NR,Nx)
 complex(kind=kind(0.d0)):: a(NR,Nstates),psi2(NR,Nx)
 complex(kind=kind(0.d0)):: sg(nr,nx), su(nr,nx)


  a = (0.d0,0.d0)
  b = 0.d0
  cx = (0.d0, 0.d0)
  ax = 0.d0

 psi2 = psi
 

  do N = 1, Nstates
   do I = 1, NR

    do J = 1, Nx
     a(i,n) = a(i,n) + ewf(j,i,n) * psi2(i,j)
    end do

   end do
  end do

  a = a * dx


  do N = 1, Nstates
   do I = 1, NR
     b(n) = b(n) + abs(a(i,n))**2
   end do
  end do

  b = b * dR


  write(500,*) sngl(time *au2fs), b(1)
  write(501,*) sngl(time *au2fs), sngl(b(2)), sngl(b(3)), sngl(b(4))


 if(mod(K,100).eq.0) then  ! writing out the probability amplitudes
   do I = 1, NR
    write(600,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,1))**2)
    write(601,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,2))**2)
   end do

   write(600,*)
   write(601,*)
 end if


   
return
end subroutine

! .......................................................................

subroutine osc_dipole(psi, d_t1, d_t2,grad)

use data_grid
use pot_param

 implicit none

 integer:: I, J
 double precision,intent(out):: d_t1, d_t2(NR)
 double precision,intent(in)::grad(nr,nx)
 complex*16,intent(in):: psi(NR,Nx)
 
 d_t1 = 0.d0
 d_t2 = 0.d0

 do i = 1, Nr
  do j = 1, Nx
   d_t2(i) = d_t2(i) + abs(psi(i,j))**2 * grad(i,j)    
  end do
   d_t2(i) = -d_t2(i) * dx
 end do 
 
 
 do i = 1, Nr
  d_t1 = d_t1 + d_t2(i)
 end do
 
  d_t1 = d_t1 * dR
  

return 
end subroutine

!........................................................................

subroutine localization(psi, time, ewf, K)

use data_grid
use data_au
use pot_param

 implicit none

 integer:: I, J, K, N

 double precision,intent(in):: time
 double precision,intent(in):: ewf(Nx,NR,Nstates)
 double precision:: B(Nstates), pl_loc(Nx,nR), neg_loc(nx,Nr)
 complex(kind=kind(0.d0)),intent(in):: psi(NR,Nx)
 complex(kind=kind(0.d0)):: psi2(NR,Nx), a(nr,Nstates)
 

 psi2 = psi
 
 pl_loc = 0.d0
 neg_loc = 0.d0
 
 do i = 1, Nr
  do j = 1, nx
   pl_loc(j,i) = 1./sqrt(2.d0)* (ewf(j,i,1) + ewf(j,i,2))
   neg_loc(j,i) = 1./sqrt(2.d0)* (ewf(j,i,1) - ewf(j,i,2))
  end do
 end do


  a = (0.d0,0.d0)
  b = 0.d0


  
   do I = 1, NR
    do J = 1, Nx
     a(i,1) = a(i,1) + pl_loc(j,i) * psi2(i,j)
     a(i,2) = a(i,2) + neg_loc(j,i) * psi2(i,j)
    end do
   end do
  

  a = a * dx

   

  do N = 1, Nstates
   do I = 1, NR
     b(n) = b(n) + abs(a(i,n))**2
   end do
  end do

  b = b * dR
  
  
 if(mod(K,100).eq.0) then  ! writing out the probability amplitudes
   do I = 1, NR
    write(504,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,1))**2)
    write(505,*) sngl(time*au2fs), sngl(R(I)*au2a),sngl(abs(a(i,2))**2)
   end do

   write(504,*)
   write(505,*)
 end if

  write(502,*) sngl(time *au2fs), b(1)
  write(503,*) sngl(time *au2fs), b(2)
 


return
end subroutine


!............... Cut off Functions ................


subroutine cofs(cof, cofR)

use data_grid
use data_au
use pot_param

 implicit none
 integer:: i, j
 double precision:: cpm, cpmR
 double precision:: cof(nx), cofr(nr)

  cof = 0.d0
  cofR = 0.d0
  cpm = 20.d0 / au2a
  cpmr = 2.d0 / au2a

!  call cutoff_cos_2d(cpm, cpmR, cof, cofR) 
 call cutoff_ex_2d(cpm,cpmR,cof,cofR)
  do I = 1, Nx/2
    cof(i) = cof(nx-i)
  end do

  do I = 1, Nx
    write(98,*) I, x(I), cof(i)
  end do



  do j = 1, NR
    write(99,*) J, R(J) * au2a, cofR(j)
  end do
 close(98, status='keep')
 close(99, status='keep')
return
end subroutine

!------------------------------------------------
subroutine cutoff_cos_2d(cpm, cpmR,cof, cofR)
use data_grid
use data_au
use pot_param

implicit none

   integer :: I,J
   double precision:: cpmR
   double precision:: cof(Nx),cpm,cofR(NR)
   do I = 1, Nx
     if(abs(x(I)).lt.abs((xend - cpm))) then
     cof(i) = 1.d0
     else
     cof(i) = cos(((x(I) - xend + cpm) / -cpm) * (0.5d0 * pi))
     cof(i) = cof(i)**2
     end if
   end do


   do j = 1, NR
    if(R(J).lt.(Rend - cpmR)) then
    cofR(j) = 1.d0
    else
    cofR(j) = dcos(((R(J) - Rend + cpmR) / -cpmR) * (0.5d0 * pi))
    cofR(j) = cofR(j)**2
    end if
   end do

 return
 end subroutine
 
 !------------------------------------------------
 subroutine cutoff_ex_2d(cpm, cpmR,cof, cofR)
 use data_grid
 use data_au
 use pot_param

  implicit none

    integer ::I, J
    double precision:: cof(Nx),cpm,c
    double precision:: cofR(NR), cpmR, cr
 open(1,file='c.out')
  c=0.300d0
  cr=3.00d0
  do I = 1, Nx
    cof(i) = 1.0d0/(1.0d0+exp(c*(x(I)-xend+cpm)))
  end do


 do j = 1, NR
   cofR(j)=1.0d0/(1.0d0+exp(cr*(R(J)-Rend+cpmR)))
   write(1,*) J, sngl(R(J)*au2a)
 end do
  write(1,*) J
  close(1)
 return
 end subroutine

