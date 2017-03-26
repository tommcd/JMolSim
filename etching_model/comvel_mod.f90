MODULE comvel_mod
   USE precision_mod
   USE global_vars
CONTAINS
  SUBROUTINE COMVEL (natom,TEMP,vxyz )
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
      real(wp):: vxyz(:,:),TEMP
      integer:: i,natom
      PI = 4.0_wp*ATAN(1.0_wp)
      PI2 = 2.0_wp*PI
      IAX1 = 987321.0_wp
      IAY1 = 199971.0_wp
      IAZ1 = 673923.0_wp
      IAX2 = 786391.0_wp
      IAY2 = 383837.0_wp
      IAZ2 = 921217.0_wp
      RMUL2 = 16807.0_wp
      RMO2 = 2147483647.0_wp
      RIMO1 = 1.0_wp/(2147483648.0_wp)
      do  i = 1,natom
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
      vxyz(i,1) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNX1))*COS(PI2*RNX2)
      vxyz(i,2) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNY1))*COS(PI2*RNY2)
      vxyz(i,3) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNZ1))*COS(PI2*RNZ2)
      end do
END SUBROUTINE
END MODULE

