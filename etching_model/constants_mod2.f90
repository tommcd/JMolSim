MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: angstrom,pi
   USE atom_types_mod
   USE keating_parameters_mod
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
!
   real(wp),parameter:: bondl_SiO = ASiO
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: AOH = bondl_OH
   real(wp),parameter:: aSiOSi = pi*(144.0_wp/180.0_wp)
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
   real(wp),parameter:: sigma(0:ntyp) = (/ sig_O_sil, &
                                          sig_Si_sil, &
                                            sigma_OH, &
                                        sig_H_OH_sil, &
                                           sig_O_H2O, &
                                           sig_H_H2O /)
   real(wp),parameter:: sigma_2(0:ntyp) = sigma*0.5_wp
   real(wp),parameter:: epsi(0:ntyp) = (/ eps_O_sil, &
                                         eps_Si_sil, &
                                       eps_O_OH_sil, &
                                       eps_H_OH_sil, &
                                          eps_O_H2O, &
                                          eps_H_H2O /)
   real(wp),parameter:: qi(0:ntyp) = (/ q_O_sil, &
                                       q_Si_sil, &
                                     q_O_OH_sil, &
                                     q_H_OH_sil, &
                                        q_O_H2O, &
                                        q_H_H2O /)
END MODULE constants_mod

