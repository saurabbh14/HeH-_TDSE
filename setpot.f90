  subroutine potential
         
use data_grid
use data_au  
use pot_param   
 
  implicit none 
  
  integer:: I, J
    
  double precision:: v12, v1e, v2e, ionicpot(Nr)
  double precision:: GM, z1, z2    


!  open(11,file='HeH+_potential_1d_0.05..0.05..102.4',status='old')
  open(10,file='gesamtpotential.out',status='unknown')
  open(23,file='ionic_pot.out',status='unknown')

  z1 = 1.d0
  z2 = 1.d0
 
 do I = 1, NR      
   v12 =(2.0d0/R(I)) + vion(I)
  do J = 1, Nx
   v1e = -1.d0 / sqrt((x(J) - mn1*R(I))**2 + en1(I))
   v2e = -1.d0 / sqrt((x(J) + mn2*R(I))**2 + en2(I))
         
  pot(I,J) = v12 + v1e + v2e 
  enddo  
 end do
 
!do J=1,Nx
!  rewind(11)
!  do I=1,NR
!  read(11,*) dummy, pot(I,J), dummy2, dummy3
!  read(11,*) dummy, dummy4, dummy2, dummy3
!  enddo
!enddo

   GM =minval(pot) 
   
   print*,'global minimum:', sngl(GM *au2eV)
   
 
  do I = 1, NR  
    do J = 1, Nx      
     if(mod(I,4).eq.0.and.mod(J,8).eq.0) then   
        write(10,*) sngl(R(I) *au2a), sngl(x(J) *au2a), sngl(pot(I,J)*au2eV)     
     end if 
    end do 
    
    if(mod(I,4).eq.0) then  
       write(10,*)      
    end if     
  end do
  
  do I = 1, NR  
    ionicpot(i) = (2) / abs(R(I))          
    write(23,*) sngl(R(I) *au2a), sngl(ionicpot(i) *au2eV)     
  end do
     
  close(10,status='keep')
  close(23,status='keep') 
  

return
end subroutine
        
            

