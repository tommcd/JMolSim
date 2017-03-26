MODULE probe_mol_init_mod
   USE precision_mod
   USE constants_mod
   USE probe_mol_mod
   implicit none
!
   type(probe_mol):: probe_SiOH4, probe_SiO4, probe_tip3p,probe_SiO4Q, probe_SiO3
   type(probe_mol):: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet1
!
CONTAINS

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types_mod
      USE constants_mod
      USE charges_mod
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigLJ_H2O*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet1,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom = (/ iOw, iHw, iHw /)
      probe_TIP3P%rad = sigLJ_2(probe_TIP3P%atom)
      probe_TIP3P%q = (/q_O_H2O, -0.5_wp*q_O_H2O, -0.5_wp*q_O_H2O /)
! Si(OH)4
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO,bondl_OH,angle_SiOH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigLJ_2(probe_SiOH4%atom)
      ! The Silicon atom should have a radius >= 1.0 A
      probe_SiOH4%rad(1) = max(1.0_wp/Angstrom,probe_SiOH4%rad(1))
      probe_SiOH4%q = qi(probe_SiOH4%atom)
      probe_SiOH4%q(1) = -4.0_wp*q_OH_sil
! SiO4, just the 4 OH Oxygens
      call new_probe_mol(probe_SiO4,4,sigLJ_OH*0.5_wp)
      call tetra_coords(bondl_SiO,probe_SiO4%r)
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigLJ_2(probe_SiO4%atom)
!
! SiO4, Si & 4 OH Oxygens
      call new_probe_mol(probe_SiO4Q,5,sigLJ_OH*0.5_wp)
      call tetra_coords5(bondl_SiO,probe_SiO4Q%r)
      probe_SiO4Q%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4Q%rad = sigLJ_2(probe_SiO4Q%atom)
      ! The Silicon atom should have a radius >= 1.0 A
      probe_SiOH4%rad(1) = max(1.0_wp/Angstrom,probe_SiOH4%rad(1))
      probe_SiO4Q%q(2:5) = q_OH_sil
      probe_SiO4Q%q(1) = -4.0_wp*q_OH_sil
! Debugging: check molecule is electroneutral
!if(sum(probe_SiO4Q%q)/=0.0_wp)then
!print *,probe_SiO4Q%q
!print *,sum(probe_SiO4Q%q)
!stop 'probe_SiO4Q%q/=0'
!end if
      allocate(probe_SiO4Q%nc(5),probe_SiO4Q%con(5,4))
      probe_SiO4Q%nc = (/ 4,1,1,1,1/)
      probe_SiO4Q%con(1,:) = (/ 2,3,4,5 /)
      probe_SiO4Q%con(2,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%con(3,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%con(4,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%con(5,:) = (/ 1,0,0,0 /)
!
! SiO3, Si & 3 OH Oxygens
      call new_probe_mol(probe_SiO3,4,sigma_OH*0.5_wp)
      probe_SiO3%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO3%rad = sigma_2(probe_SiO3%atom)
      allocate(probe_SiO3%nc(4),probe_SiO3%con(4,4))
      probe_SiO3%nc = (/ 3,1,1,1 /)
      probe_SiO3%con(1,:) = (/ 2,3,4,0 /)
      probe_SiO3%con(2,:) = (/ 1,0,0,0 /)
      probe_SiO3%con(3,:) = (/ 1,0,0,0 /)
      probe_SiO3%con(4,:) = (/ 1,0,0,0 /)
      probe_SiO3%r(:,1) = probe_SiO4Q%r(:,1)
      probe_SiO3%r(:,2) = probe_SiO4Q%r(:,3)
      probe_SiO3%r(:,3) = probe_SiO4Q%r(:,4)
      probe_SiO3%r(:,4) = probe_SiO4Q%r(:,5)

! O2
      call new_probe_mol(probe_O2,2,sig_O_O2*0.5_wp)
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_O2*0.5_wp /)
! N2
      call new_probe_mol(probe_N2,2,sig_N_N2*0.5_wp)
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_N2*0.5_wp /)
! CO2
      call new_probe_mol(probe_CO2,3,sig_O_CO2*0.5_wp)
      probe_CO2%rad(1) = sig_C_CO2*0.5_wp
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp,-bondl_CO2 /)
!
! Write out details of the probe molecules
! Si(OH)4
      write(*,*) 'SiOH4'
      write(*,*) probe_SiOH4%n
      write(*,'(10i6)') probe_SiOH4%atom
      write(*,'(9f12.6)') probe_SiOH4%rad
      write(*,'(9f12.6)') probe_SiOH4%q
      write(*,'(3f12.6)') probe_SiOH4%r
! H2O
      write(*,*) 'H2O TIP3P'
      write(*,*) probe_TIP3P%n
      write(*,'(10i6)') probe_TIP3P%atom
      write(*,'(9f12.6)') probe_TIP3P%rad
      write(*,'(9f12.6)') probe_TIP3P%q
      write(*,'(3f12.6)') probe_TIP3P%r
! SiO4
      write(*,*) 'SiO4 ... just OH'
      write(*,*) probe_SiO4%n
      write(*,'(9f12.6)') probe_SiO4%rad
      write(*,'(3f12.6)') probe_SiO4%r
! Si(OH)4
      write(*,*) 'SiO4'
      write(*,*) probe_SiO4Q%n
      write(*,'(9i6)') probe_SiO4Q%atom
      write(*,'(9f12.6)') probe_SiO4Q%rad
      write(*,'(9f12.6)') probe_SiO4Q%q
      write(*,'(3f12.6)') probe_SiO4Q%r
! Si(OH)3
      write(*,*) 'SiO3 ... just OH'
      write(*,*) probe_SiO3%n
      write(*,'(9i6)') probe_SiO3%atom
      write(*,'(9f12.6)') probe_SiO3%rad
      write(*,'(3f12.6)') probe_SiO3%r

! O2
      write(*,*) 'O2'
      write(*,*) probe_O2%n
      write(*,'(9f12.6)') probe_O2%rad
      write(*,'(3f12.6)') probe_O2%r
! N2
      write(*,*) 'N2'
      write(*,*) probe_N2%n
      write(*,'(9f12.6)') probe_N2%rad
      write(*,'(3f12.6)') probe_N2%r
! CO2
      write(*,*) 'CO2'
      write(*,*) probe_CO2%n
      write(*,'(9f12.6)') probe_CO2%rad
      write(*,'(3f12.6)') probe_CO2%r
   END SUBROUTINE init_probe_mols


   SUBROUTINE h2o_coords(bondl,ang,r)
      real(wp),intent(in):: bondl,ang
      real(wp),intent(out):: r(3,3)
! Coordinates for water molecule
      r(1:3,1) = 0.0_wp
      r(1:3,2) = (/ 0.0_wp, 0.0_wp, bondl /)
      r(1:3,3) = (/ 0.0_wp, bondl*sin(ang), bondl*cos(ang) /)
   END SUBROUTINE

END MODULE probe_mol_init_mod

