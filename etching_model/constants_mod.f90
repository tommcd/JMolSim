MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars, only: angstrom,pi
   USE keating_parameters
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
!
   real(wp),parameter:: AOH = bondl_OH
!  Dreiding
   real(wp),parameter:: KOH = 3035.2_wp    ! ((700 kcal/mol)/A2)   eV nm^-2
   real(wp),parameter:: KSiOH = 4.336_wp   ! ((100 kcal/mol)/rad2)
!
! the following are now taken from Seaton
!  real(wp),parameter:: sigma_Obridge,sigma_OH,sigma_O,sigma_Si,sigma_C,sigma_N
!  real(wp),parameter:: bondl_O2,bondl_N2,bondl_SiO,bondl_CO2
!  real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
! united atom SiOH4
   real(wp),parameter:: r_tet = bondl_SiO + sigma_OH*0.5_wp
!
   real(wp),parameter:: sigma(0:3) = (/ sig_O_sil, &
                                        sig_Si_sil,&
                                        sigma_OH,  &
                                        sig_H_OH_sil /)
   real(wp),parameter:: sigma_2(0:3) = sigma*0.5_wp
   real(wp),parameter:: epsi(0:3) = (/ eps_O_sil, &
                                       eps_Si_sil, &
                                       eps_O_OH_sil, &
                                       eps_H_OH_sil /)
END MODULE constants_mod

