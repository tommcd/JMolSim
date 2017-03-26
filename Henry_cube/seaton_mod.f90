
MODULE seaton_mod
   USE precision_mod
   USE global_vars, only: angstrom, K_ev, qstar, pi
   implicit none
! atom types
   integer,parameter:: iO_Sil=0
   integer,parameter:: iSi_Sil=1
   integer,parameter:: iO_OH_sil=2
   integer,parameter:: iH_OH_sil=3
   integer,parameter:: iO_H2O=4
   integer,parameter:: iH_H2O=5
   integer,parameter:: iO_CO2=6
   integer,parameter:: iC_CO2=7
   integer,parameter:: iN_N2=8
   integer,parameter:: iO_O2=9
   integer,parameter:: iN_N2charge=10
   integer,parameter:: iO_O2charge=11
   integer,parameter:: ntyplj = 11  ! number of LJ types
   real(wp),parameter:: eps_O_sil   =185.0_wp*K_ev, sig_O_sil  = 2.708_wp/angstrom,q_O_sil = -0.64025_wp*qstar
   real(wp),parameter:: eps_Si_sil  =  0.0_wp*K_ev, sig_Si_sil = 0.0_wp/angstrom,  q_Si_sil= -2.0_wp*q_O_sil
   real(wp),parameter:: eps_O_OH_sil=185.0_wp*K_ev, sig_O_OH_sil= 3.0_wp/angstrom, q_O_OH_sil= -0.533_wp*qstar
   real(wp),parameter:: eps_H_OH_sil=  0.0_wp*K_ev, sig_H_OH_sil= 0.0_wp/angstrom, q_H_OH_sil=  0.206_wp*qstar
   real(wp),parameter:: eps_H_H2O   =  0.0_wp*K_ev,  sig_H_H2O = 0.0_wp/angstrom,   q_H_H2O=  0.417_wp*qstar
   real(wp),parameter:: eps_O_H2O   = 76.58_wp*K_ev, sig_O_H2O = 3.1506_wp/angstrom,q_O_H2O= -2.0_wp*q_H_H2O
   real(wp),parameter:: eps_O_CO2   = 82.997_wp*K_ev,sig_O_CO2 = 3.064_wp/angstrom, q_O_CO2= -0.33225_wp*qstar
   real(wp),parameter:: eps_C_CO2   = 29.999_wp*K_ev,sig_C_CO2 = 2.785_wp/angstrom, q_C_CO2= -2.0_wp*q_O_CO2
   real(wp),parameter:: eps_N_N2    = 34.897_wp*K_ev, sig_N_N2 = 3.3211_wp/angstrom, q_N_N2= -0.5475_wp*qstar
   real(wp),parameter:: eps_O_O2    = 43.183_wp*K_ev, sig_O_O2 = 3.1062_wp/angstrom, q_O_O2= -0.3577_wp*qstar
!  bond lengths
   real(wp),parameter:: bondl_SiO = 1.60_wp/angstrom
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
! LJ parameters
   real(wp),parameter:: epsi(0:ntyplj)=(/ &
                        eps_O_sil, &
                        eps_Si_sil, &
                        eps_O_OH_sil, &
                        eps_H_OH_sil, &
                        eps_O_H2O, &
                        eps_H_H2O, &
                        eps_O_CO2, &
                        eps_C_CO2, &
                        eps_N_N2, &
                        eps_O_O2, &
                        0.0_wp, &
                        0.0_wp /)
   real(wp),parameter:: sigi(0:ntyplj)=(/ &
                        sig_O_sil, &
                        sig_Si_sil, &
                        sig_O_OH_sil, &
                        sig_H_OH_sil, &
                        sig_O_H2O, &
                        sig_H_H2O, &
                        sig_O_CO2, &
                        sig_C_CO2, &
                        sig_N_N2, &
                        sig_O_O2, &
                        0.0_wp, &
                        0.0_wp /)
   real(wp),parameter:: sigi2(0:ntyplj)=sigi*0.5_wp
! Charges
   real(wp),parameter:: q_seaton(0:11)=(/ &
                        q_O_sil, &
                        q_Si_sil, &
                        q_O_OH_sil, &
                        q_H_OH_sil, &
                        q_O_H2O, &
                        q_H_H2O, &
                        q_O_CO2, &
                        q_C_CO2, &
                        q_N_N2, &
                        q_O_O2, &
                -2.0_wp*q_N_N2, &
                -2.0_wp*q_O_O2 /)
END MODULE seaton_mod

