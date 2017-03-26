
MODULE global_vars
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: length_m = angstrom*1.0e-10_wp
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: amu = 1.66053873e-27_wp
   real(wp),parameter:: eV = 1.602176462e-19_wp
   real(wp),parameter:: k_boltz = 1.3806503e-23_wp
   real(wp),parameter:: K_ev = k_boltz/eV
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: time_s,dt_fs
   real(wp):: del_rxn, e_activ, etot
   integer:: natom,nseed,natom_max
   integer:: nattached, nrelax
END MODULE

