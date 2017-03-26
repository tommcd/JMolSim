
module Lj_el_mod
    use precision_mod
    use coordinates_mod 
!    use Keating_parameters
    use charges_mod
    use atom_types
    USE seaton_mod
    USE constants_mod
    implicit none
    save
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: vdwcut(:,:),vdwsf(:,:)
    real(wp),parameter:: rcut = 12.0_wp/angstrom
    real(wp),parameter:: rcut2 = rcut**2
    real(wp),parameter:: theta_CO2 =3.1415926535 
    real(wp),parameter:: K_CO2 = 5.331    
    real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3,rl1(3),rl1sq,rl2(3),rl2sq
    real(wp):: theta,dist,rr13,df(3),rcut6
    

 
CONTAINS
     

   SUBROUTINE LJ_INIT

      real(wp):: eij,rstij
      integer:: i,j
    allocate( aij(0:ntyplj,0:ntyplj) )
    allocate( bij(0:ntyplj,0:ntyplj) )
print *,'ntyplj = ',ntyplj
    allocate( vdwcut(0:ntyplj,0:ntyplj) )
    allocate( vdwsf(0:ntyplj,0:ntyplj) ) 

      do j = 0,ntyplj    ! lj aij,bij parameters & shifted terms
      do i = 0,ntyplj

        rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
        eij = sqrt(epsi(i)*epsi(j))
        rstij = (sigi(i)+sigi(j))*0.5
        aij(i,j) = 4.0*(rstij**12)*eij
        bij(i,j) = 4.0*(rstij**6)*eij
        vdwcut(i,j) = rcut6*(aij(i,j)*rcut6-bij(i,j)) 
        vdwsf(i,j) = rcut6*(-12.0_wp*aij(i,j)*rcut6+6.0*bij(i,j))/rcut
!        write(*,*)i,j,aij(i,j),bij(i,j),vdwcut(i,j),vdwsf(i,j)
      end do
      end do

   END SUBROUTINE LJ_INIT

   SUBROUTINE FORCE_ENERGY_LJ_EL_O2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      REAL(wp),intent(out):: ULJEL
      real(wp):: rr(3),r2,r1,rr6,dele,rr7,rcut7,rcut13
      integer:: L,M,n,i
      ULJEL=0.0_wp
      outer_atom_loop: DO L = atomfirst,atomlast
      atom_loop: DO M = sysfirst,syslast 
         Imige_Ox_xyz(M,1,:) = Ox_xyz(M,1,:)
         rr(1:3) = Ox_xyz(M,2,:) - Ox_xyz(M,1,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         Imige_Ox_xyz(M,2,:) = Ox_xyz(M,1,:) + rr(:) 
         Imige_Ox_xyz(M,3,:) = (Imige_Ox_xyz(M,1,:) + Imige_Ox_xyz(M,2,:))/2.0_wp

       DO n = 1,2 
         rr(1:3) = Ox_xyz(M,n,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)       
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         IF (r2 < rcut2) THEN
            rr6 = 1.0_wp/(r2*r2*r2)
            rr13 = rr6*rr6/r1
            rr7 = rr6/r1
            rcut7 = 1.0_wp/((rcut2*rcut2*rcut2)*rcut)
            rcut13 = rcut7*rcut7*rcut
            dele = rr6*(aij(atom(L),Ox_atom(M,n))*rr6-bij(atom(L),Ox_atom(M,n))) &
                 - vdwcut(atom(L),Ox_atom(M,n)) - vdwsf(atom(L),Ox_atom(M,n))*(r1-rcut)
            df = (6.0_wp*bij(atom(L),Ox_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),Ox_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            Ox_fxyz(M,n,1:3) = Ox_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + dele
         END IF

         rr(1:3) = Imige_Ox_xyz(M,n,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
            df = - charge(L)*Ox_gas_charge(M,n)*rr(1:3)/(r2*r1) 
            fxyz(L,1:3) = fxyz(L,1:3) + df
            Ox_fxyz(M,n,1:3) = Ox_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,n)/r1
!           
      END DO 
!               
         rr(1:3) = Imige_Ox_xyz(M,3,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         df = - charge(L)*Ox_gas_charge(M,3)*rr(1:3)/(r2*r1) 
            fxyz(L,1:3) = fxyz(L,1:3) + df
            Ox_fxyz(M,1,1:3) = Ox_fxyz(M,1,1:3) - 0.5_wp*df 
            Ox_fxyz(M,2,1:3) = Ox_fxyz(M,2,1:3) - 0.5_wp*df
         ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,3)/r1
        
      
      END DO atom_loop
      END DO outer_atom_loop
      energy = energy + ULJEL 
   END SUBROUTINE
   
   SUBROUTINE FORCE_ENERGY_LJ_EL_N2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      REAL(wp),intent(out):: ULJEL
      real(wp):: rr(3),r2,r1,rr6,dele,rr7,rcut7,rcut13
      integer:: L,M,i,n
      ULJEL=0.0_wp
      outer_atom_loop: DO L = atomfirst,atomlast
      atom_loop: DO M = sysfirst,syslast 
         Imige_N2_xyz(M,1,:) = N2_xyz(M,1,:)
         rr(1:3) = N2_xyz(M,2,:) - N2_xyz(M,1,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         Imige_N2_xyz(M,2,:) = N2_xyz(M,1,:) + rr(:)
         Imige_N2_xyz(M,3,:) = (Imige_N2_xyz(M,1,:) + Imige_N2_xyz(M,2,:))/2.0_wp

       DO n = 1,2  
         rr(1:3) = N2_xyz(M,n,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)                
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         IF (r2 < rcut2) THEN
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            rcut7 = 1.0_wp/(((rcut**3)**2)*rcut)
            rcut13 = rcut7*rcut7*rcut
            dele = rr6*(aij(atom(L),N2_atom(M,n))*rr6-bij(atom(L),N2_atom(M,n))) &
                 - vdwcut(atom(L),N2_atom(M,n)) - vdwsf(atom(L),N2_atom(M,n))*(r1-rcut)
            df = (6.0_wp*bij(atom(L),N2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),N2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            N2_fxyz(M,n,1:3) = N2_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + dele
         END IF

         rr(1:3) = Imige_N2_xyz(M,n,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
            df = - charge(L)*N2_gas_charge(M,n)*rr(1:3)/(r2*r1) 
            fxyz(L,1:3) = fxyz(L,1:3) + df
            N2_fxyz(M,n,1:3) = N2_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + charge(L)*N2_gas_charge(M,n)/r1
           
      END DO       
         rr(1:3) = Imige_N2_xyz(M,3,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         df = - charge(L)*N2_gas_charge(M,3)*rr(1:3)/(r2*r1) 
            fxyz(L,1:3) = fxyz(L,1:3) + df
            N2_fxyz(M,1,1:3) = N2_fxyz(M,1,1:3) - 0.5_wp*df 
            N2_fxyz(M,2,1:3) = N2_fxyz(M,2,1:3) - 0.5_wp*df

         ULJEL = ULJEL + charge(L)*N2_gas_charge(M,3)/r1
        
      
      END DO atom_loop
      END DO outer_atom_loop
      energy = energy + ULJEL 
   END SUBROUTINE

   subroutine O2_streching 
!     Force routine for Keating bond streching
      integer:: i,n
      do i = 1,n_O2
        rl1 = Ox_xyz(i,1,:)-Ox_xyz(i,2,:)
        !call pbc(rl1)
        rl1 = rl1-boxl*anint(rl1*boxli)

        rl1sq = sqrt(dot_product(rl1,rl1)) 

        Ox_fxyz(i,1,1:3) = Ox_fxyz(i,1,1:3)-KOO*(rl1sq-AOO)*rl1(1:3)/rl1sq
        Ox_fxyz(i,2,1:3) = Ox_fxyz(i,2,1:3)+KOO*(rl1sq-AOO)*rl1(1:3)/rl1sq

        energy = energy + 0.5_wp*KOO*(rl1sq-AOO)**2
        
      end do
    end subroutine O2_streching

   subroutine N2_streching 
!     Force routine for Keating bond streching
      integer:: i,n
      do i = 1,n_N2
        rl1 = N2_xyz(i,1,:)-N2_xyz(i,2,:)
        !call pbc(rl1)
        rl1 = rl1-boxl*anint(rl1*boxli)

        rl1sq = sqrt(dot_product(rl1,rl1)) 

        N2_fxyz(i,1,1:3) = N2_fxyz(i,1,1:3)-KNN*(rl1sq-ANN)*rl1(1:3)/rl1sq
        N2_fxyz(i,2,1:3) = N2_fxyz(i,2,1:3)+KNN*(rl1sq-ANN)*rl1(1:3)/rl1sq

        energy = energy + 0.5_wp*KNN*(rl1sq-ANN)**2
        
      end do
    end subroutine

    SUBROUTINE FORCE_ENERGY_CO2_LJ_EL(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      REAL(wp),intent(out):: ULJEL
      real(wp):: rr(3),r2,r1,rr6,dele,rr7,rcut7,rcut13
      integer:: L,M,i,n
      ULJEL=0.0_wp
      outer_atom_loop: DO L = atomfirst,atomlast
      atom_loop: DO M = sysfirst,syslast 
         Imige_CO2_xyz(M,1,:) = CO2_xyz(M,1,:)
         rr(1:3) = CO2_xyz(M,2,:) - CO2_xyz(M,1,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         Imige_CO2_xyz(M,2,:) = CO2_xyz(M,1,:) + rr(:)
         rr(1:3) = CO2_xyz(M,3,:) - CO2_xyz(M,1,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
         Imige_CO2_xyz(M,3,:) = CO2_xyz(M,1,:) + rr(:)         

         DO n = 1,3 
         rr(1:3) = CO2_xyz(M,n,:) - rxyz(L,:)
         rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)                
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         IF (r2 < rcut2) THEN
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            rcut7 = 1.0_wp/(((rcut**3)**2)*rcut)
            rcut13 = rcut7*rcut7*rcut
            dele = rr6*(aij(atom(L),CO2_atom(M,n))*rr6-bij(atom(L),CO2_atom(M,n))) &
                 - vdwcut(atom(L),CO2_atom(M,n)) - vdwsf(atom(L),CO2_atom(M,n))*(r1-rcut)
            df = (6.0_wp*bij(atom(L),CO2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),CO2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            CO2_fxyz(M,n,1:3) = CO2_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + dele
         END IF
            rr(1:3) = Imige_CO2_xyz(M,n,:) - rxyz(L,:)
            rr(1:3) = rr(1:3) - boxl*ANINT(rr(1:3)*boxli)
            r2 = rr(1)**2+rr(2)**2+rr(3)**2
            r1 = sqrt(r2)

            df = - charge(L)*CO2_gas_charge(M,n)*rr(1:3)/(r2*r1) 
            fxyz(L,1:3) = fxyz(L,1:3) + df
            CO2_fxyz(M,n,1:3) = CO2_fxyz(M,n,1:3) - df 
            ULJEL = ULJEL + charge(L)*CO2_gas_charge(M,n)/r1

         END DO 
                  
      END DO atom_loop
      END DO outer_atom_loop
      energy = energy + ULJEL 
   END SUBROUTINE

   subroutine CO2_streching 
!     Force routine for Keating bond streching
      integer:: i,n
      do i = 1,n_CO2
        dist = ACO2
        rl1 = CO2_xyz(i,2,:)-CO2_xyz(i,1,:)
        !call pbc(rl1)
        rl1 = rl1-boxl*anint(rl1*boxli)

        rl1sq = sqrt(dot_product(rl1,rl1)) 

        CO2_fxyz(i,2,1:3) = CO2_fxyz(i,2,1:3)-KCO2*(rl1sq-dist)*rl1(1:3)/rl1sq
        CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3)+KCO2*(rl1sq-dist)*rl1(1:3)/rl1sq

        energy = energy + 0.5_wp*KCO2*(rl1sq-dist)**2
        
        rl1 = CO2_xyz(i,3,:)-CO2_xyz(i,1,:)
        rl1 = rl1-boxl*anint(rl1*boxli)
        rl1sq = sqrt(dot_product(rl1,rl1))
        CO2_fxyz(i,3,1:3) = CO2_fxyz(i,3,1:3)-KCO2*(rl1sq-dist)*rl1(1:3)/rl1sq
        CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3)+KCO2*(rl1sq-dist)*rl1(1:3)/rl1sq

        energy = energy + 0.5_wp*KCO2*(rl1sq-dist)**2        
      end do
    end subroutine

    subroutine CO2_angle_bending
!     Force routine for CO2 bond angle bending
      integer:: i,n
      do i = 1,n_CO2 
        rl1 =CO2_xyz(i,2,1:3)-CO2_xyz(i,1,1:3)
        rl2 =CO2_xyz(i,3,1:3)-CO2_xyz(i,1,1:3)

!        !call pbc(rl1)
!        !call pbc(rl2)
        rl1 = rl1-boxl*anint(rl1*boxli)
        rl2 = rl2-boxl*anint(rl2*boxli)
!
        rl1sq = sqrt(dot_product(rl1,rl1)) 
        rl2sq = sqrt(dot_product(rl2,rl2))
        rl1l2 = dot_product(rl1,rl2)
        cos_theta = (rl1l2)/(rl1sq*rl2sq)
        theta = acos(cos_theta)

       IF (sqrt(1-cos_theta**2) == 0.0_wp) then 
        CO2_fxyz = CO2_fxyz
         else
        coef = K_CO2*(theta-theta_CO2)*(-1.0_wp/sqrt(1-cos_theta**2))
        coef1 = rl1sq**3*rl2sq
        coef2 = rl1sq*rl2sq 
        coef3 = rl1sq*rl2sq**3 



         CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3)-coef*(-(rl1+rl2)/coef2 &
                                        + rl1l2*rl1/coef1 &
                                        + rl1l2*rl2/coef3)
         CO2_fxyz(i,2,1:3) = CO2_fxyz(i,2,1:3)-coef*(rl2/coef2-rl1l2*rl1/coef1)
         CO2_fxyz(i,3,1:3) = CO2_fxyz(i,3,1:3)-coef*(rl1/coef2-rl1l2*rl2/coef3)
         
              
         END IF   
          energy = energy + 0.5_wp* K_CO2*(theta-theta_CO2)**2
        
      end do

   end subroutine    
END MODULE Lj_el_mod

! electrostatic
!          chargerij = charge(i)*charge(j)*rabi
!          emel = emel+chargerij








