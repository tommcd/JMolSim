module comvel_mod
   use precision_mod
   use global_vars
contains
  SUBROUTINE COMVEL (natom,n_CO2,n_O2,n_N2,TEMP,vxyz,Ox_vxyz,N2_vxyz,CO2_vxyz)    
!   *******************************************************************
!   ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
!   **                                                               **
!   ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (UNIT) MASS.**
!   ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
!   ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
!   **                                                               **
!   ** ROUTINE REFERENCED:                                           ** 
!   **                                                               **
!   ** REAL FUNCTION GAUSS ( DUMMY )                                 **
!   **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
!   **    DISTRIBUTION WITH ZERO MEAN AND UNIT VARIANCE.             **
!   *******************************************************************
!
! INITIAL VELOCITY FOR MOVING PARTICLE        
!
      real(wp) PI,PI2,RMUL2,IAX1,IAX2,RMO2,IAY1,IAY2,IAZ1,IAZ2
      real(wp):: RIMO1,RNX1,RNX2,RNY1,RNY2,RNZ1,RNZ2
      real(wp):: vxyz(:,:),Ox_vxyz(:,:,:),N2_vxyz(:,:,:),CO2_vxyz(:,:,:),TEMP
      integer:: i,natom,n_CO2,n_Ox,n_N2
      PI = 4.0D0*DATAN(1.0D0)
      PI2 = 2.0*PI
      IAX1=987321.0d0
      IAY1=199971.0d0
      IAZ1=673923.0d0
      IAX2=786391.0d0
      IAY2=383837.0d0
      IAZ2=921217.0d0
      RMUL2=16807.0d0
      RMO2=2147483647.0d0
      RIMO1=1.0d0/(2147483648.0d0) 
      do  i=1,natom + 2*n_O2 + 2*n_N2 + 3*n_CO2                                                             
! INITIAL VELOCITY FROM MAXWELLIAN DISTRIBUTION                                                          
      IAX1 = MOD(RMUL2*IAX1,RMO2)                      
      IAX1 = MOD(RMUL2*IAX1,RMO2)                
      IAX2 = MOD(RMUL2*IAX2,RMO2)                        
      IAX2 = MOD(RMUL2*IAX2,RMO2)            
      IAY1 = MOD(RMUL2*IAY1,RMO2)            
      IAY1 = MOD(RMUL2*IAY1,RMO2)         
      IAY2 = MOD(RMUL2*IAY2,RMO2)                        
      IAY2 = MOD(RMUL2*IAY2,RMO2)                   
      IAZ1 = MOD(RMUL2*IAZ1,RMO2)                    
      IAZ1 = MOD(RMUL2*IAZ1,RMO2)          
      IAZ2 = MOD(RMUL2*IAZ2,RMO2)                         
      IAZ2 = MOD(RMUL2*IAZ2,RMO2)                                                                            
      RNX1 = IAX1*RIMO1                                       
      RNX2 = IAX2*RIMO1                                           
      RNY1 = IAY1*RIMO1                                             
      RNY2 = IAY2*RIMO1                              
      RNZ1 = IAZ1*RIMO1                                   
      RNZ2 = IAZ2*RIMO1                      
      vxyz(i,1)=dsqrt(TEMP*K_ev)*DSQRT(-2.D0*DLOG(RNX1))*DCOS(PI2*RNX2)                
      vxyz(i,2)=dsqrt(TEMP*K_ev)*DSQRT(-2.D0*DLOG(RNY1))*DCOS(PI2*RNY2)                  
      vxyz(i,3)=dsqrt(TEMP*K_ev)*DSQRT(-2.D0*DLOG(RNZ1))*DCOS(PI2*RNZ2)                  
      end do
         j = 1
      do i=1,n_O2
       Ox_vxyz(i,1,:) = vxyz(natom+j,:)
       Ox_vxyz(i,2,:) = vxyz(natom+j+1,:)
        j = j + 2
      end do
         j = 1
      do i=1,n_N2
       N2_vxyz(i,1,:) = vxyz(natom+2*n_O2+j,:)
       N2_vxyz(i,2,:) = vxyz(natom+2*n_O2+j+1,:)
        j = j + 2
      end do
         j = 1
      do i=1,n_CO2
       CO2_vxyz(i,1,:) = vxyz(natom+2*n_O2+2*n_N2+j,:)
       CO2_vxyz(i,2,:) = vxyz(natom+2*n_O2+2*n_N2+j+1,:)
       CO2_vxyz(i,3,:) = vxyz(natom+2*n_O2+2*n_N2+j+2,:)       
        j = j + 3
      end do
      
end subroutine
end module

