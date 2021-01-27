module data_grid
 integer:: Nt, Nstates, Vstates
 integer, parameter:: Nx=1024,NR=512 !changed back  !This program has been edited to check momentum grid dependance on number of grid point
 !reverted  back the changes while running calculations(Change date: 11/03/20)
 double precision:: RI, kappa,dR, R(NR), x(Nx)
 double precision:: vion(NR), en1(NR), en2(NR)
 double precision:: dx, dt, xeq
 double precision:: dpr, dpx, fwhm, t_start
 double precision:: e01, e02, phi1, phi2
 double precision:: omega1, omega2
 double precision:: Px(Nx),PR(NR)
 double precision:: Pot(NR,Nx)
 double precision:: mn, mn1, mn2, m1, m2 !all relevent mass veriables
end module

module data_au
 double precision,parameter:: au2a=0.52917706d0  ! length in a.u. to Angstrom
 double precision,parameter:: cm2au=4.5554927d-6 ! energy in wavenumbers to a.u.
 double precision,parameter:: au2fs=0.024        ! time in a.u. to femtoseconds
 double precision,parameter:: j2eV = 6.242D18    ! transforms energy in J to eV
 double precision,parameter:: au2eV=27.2116d0    ! energy in a.u. to eV
 double precision,parameter:: i2au=2.0997496D-9
 double precision,parameter:: e02au=5.142206707e11
 double precision,parameter:: pi=3.141592653589793d0       ! that's just pi
 double precision,parameter:: mass=1836.15d0    ! reduced mass of deuterium
 double precision,parameter:: me =1.d0          ! mass of the electron
 double precision:: m_eff, m_red                ! effecitve mass of el. and nucl.
 complex*16,parameter:: im=(0.d0,1.d0) 
end module

module pot_param
use data_au
 double precision:: R0     ! Grid-Parameter, start..
 double precision,parameter:: x0 = -100.d0!/au2a
 double precision::Rend   !..and end
 double precision,parameter:: xend = 100.d0!/au2a
end module pot_param

module FFTW3
  use, intrinsic :: iso_c_binding
  ! include '/usr/include/fftw3.f03'                                        ! Desktop packet
  ! include '/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/fftw3.f03' ! ARA cluster
  include '/usr/include/fftw3.f03'
end module



program metiu

use data_grid
 implicit none
 integer:: I, J
 
 Real*4 st, ft
 double precision, allocatable, dimension(:,:,:):: mu_all
 double precision, allocatable, dimension(:,:):: adb
 double precision, allocatable, dimension(:,:,:):: ewf
 double precision, allocatable, dimension(:,:):: chi0
  call cpu_time(st)
  call input
  call p_grid
  
 allocate(adb(NR,Nstates),ewf(Nx,NR,Nstates),mu_all(NR,Nstates,Nstates))
 allocate(chi0(nr,vstates))
  

  call potential
    
  ewf = 0.d0
  adb = 0.d0 

  call adiabatic_surface(adb, ewf, mu_all)  

  call nuclear_wavefkt(adb,chi0)
!  call propagation_1D(adb, mu_all, chi0)
  deallocate (adb, mu_all)

 call propagation_2D(ewf,chi0)          
    
 deallocate (ewf, chi0)

 call cpu_time(ft)
 print*,'Run time=', ft-st, 'seconds' 
 
stop  
end program


! _______________ Subroutines __________________________________________________


subroutine input

use data_grid
use pot_param
implicit none
double precision:: tp, lambda1, lambda2
 real*4 :: dummy, dummy2, dummy3, dummy4
 integer:: i, j

 open(11,file='HeH+_potential_1d_0.05..0.05..102.4',status='old')
 do i=1,NR
 read(11,*) R(I), vion(I), en1(I), en2(I)
 read(11,*) dummy, dummy2, dummy3, dummy4
 enddo
 close(11)

 open(10,file='input',status='old')
 
  read(10,*) m1, m2               ! m1=mass1(helium), m2=mass2(hydrogen)
  read(10,*) dt                   ! dt = time step in a.u.
  read(10,*) Nt                   ! Nt = number of time steps.	
  read(10,*) Nstates              ! Nstates = Number of calculated excited states.
  read(10,*) Vstates              ! Vstates = Number of calculated vibrational states.
  read(10,*) RI                   ! RI = center of initial Gaussian in 10^-10 m.
  read(10,*) kappa                ! kappa = width of initial Gaussian 
  read(10,*) lambda1, lambda2     ! wavelength of pulse
  read(10,*) tp, t_start, E01, E02
  read(10,*) phi1, phi2                  ! carrier-envelope phase
  
  R0=R(1)
  Rend=R(NR)
  dR =R(2)-R(1)! (Rend - R0) / (NR - 1)
  dx = (xend - x0) / (Nx - 1)
  dpx = (2.d0 * pi) / (dx * Nx)      
  dpr = (2.d0 * pi) / (dR * NR) 


 !Grids: co-ordinate space

  do J=1,Nx
   x(J) = x0 + (J - 1) * dx
  enddo
!--------------------------
!Masses:
  m1=m1*mass
  m2=m2*mass 
  mn=m1+m2
  m_red = m1*m2/(m1+m2)
!  m_eff= (m1*me+m2*me+m1*m2)/(m1*m2*me)
  m_eff = (m1+m2) / (m1+m2+1.0d0)!(4.d0 * me * mass) / (2.d0 * mass + me) 
  mn1=m1/mn
  mn2=m2/mn
!----------------------------  
  RI= RI / au2a   
  
  tp = tp / au2fs  
  t_start = t_start / au2fs   
  fwhm = (4.d0 * dlog(2.d0)) / tp**2 
  omega1=(1.d0 / (lambda1 * 1.d-7)) *cm2au
  omega2=(1.0d0/ (lambda2 *1.d-7))* cm2au
  phi1=phi1*pi
  phi2=phi2*pi
  print*,'_________________________'
  print*
  print*,'Parameters'
  print*,'_________________________'
  print*
  print*,'dt = ', SNGL(dt), 'a.u.'
  print*,'dx = ', SNGL(dx), 'a.u.'
  print*,'dR = ', SNGL(dR), 'a.u.'
  print*,'dpx = ', SNGL(dpx), 'a.u.'
  print*,'dPR = ', SNGL(dpR), 'a.u.'
  print*,'RI=', sngl(RI), 'a.u.'
  print*,'R0=', sngl(R0), 'a.u.', 'Rend=',sngl(Rend), 'a.u.'
  print*,'Wavelength 1 =', sngl(lambda1), 'nm'
  print*,'Phase 1 =', sngl(phi1)
  print*,'Field strength =', sngl(e01), 'a.u.', sngl(e01*e02au), 'V/m'
  print*,'Intensity =', sngl(e01**2*3.509e16), 'W/cm2'
  print*,'Wavelength 2 =', sngl(lambda2), 'nm'
  print*,'Phase 2 =', sngl(phi2)
  print*,'Field strength =', sngl(e02), 'a.u.', sngl(e02*e02au), 'V/m'
  print*,'Intensity =', sngl(e02**2*3.509e16), 'W/cm2'
  print*,'Wave duration =', sngl(tp*au2fs), 'fs'
  print*
  print*,'__________________________'
  print*   
                 
 end subroutine
 
!...................... Impulsgrid......................
 
subroutine p_grid

use data_grid
use pot_param
use data_au
 implicit none 
 integer:: I 
  
      
  do I = 1, Nx  
    if (I.le.(Nx / 2)) then    
    Px(I) = (I - 1) * dpx    
    else    
    Px(I) = - (Nx + 1 - I) * dpx    
    end if    
  end do
    print*, 'x0=', x0, 'xend=', xend
    print*, 'Px0=', Px((Nx/2)+1), 'Pxend=', Px(Nx/2) 
  
  do I = 1, NR  
    if (I.le.(NR / 2)) then    
    PR(I) = (I - 1) * dpR    
    else    
    PR(I) = - (NR + 1 - I) * dpR    
    end if
  end do
    print*, 'R0=', R(1), 'Rend=', R(NR)
    print*, 'PR0=', PR((NR/2)+1), 'PRend=', PR(NR/2)
        

  return
end subroutine  

