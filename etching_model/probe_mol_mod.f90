
MODULE probe_mol_mod
   USE precision_mod
   USE constants_mod
   implicit none
!
!  definition of probe_mol derived type
!  the allocatable component requires a standard extension to F95
!  described in the technical report ISO/IEC TR 15581
!
   TYPE probe_mol
      integer:: n
      real(wp),allocatable:: r(:,:)
      integer,allocatable:: atom(:)
      real(wp),allocatable:: rad(:),q(:)
   END TYPE probe_mol
!
   type(probe_mol):: probe_SiOH4, probe_SiO4, probe_tip3p
   type(probe_mol):: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet
!
CONTAINS

   SUBROUTINE new_probe_mol(P,n,radius,atomtype)
      type(probe_mol),intent(out):: P
      integer,intent(in):: n
      real(wp),intent(in),optional:: radius
      integer,intent(in),optional:: atomtype
      P%n = n
      allocate(P%r(3,n),P%atom(n),P%rad(n),P%q(n))
      P%r = 0.0_wp
      if (present(radius)) P%rad = radius
      if (present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types
      USE seaton_mod
      USE charges_mod
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigma_H2O*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom = (/ iO_H2O, iH_H2O, iH_H2O /)
      probe_TIP3P%rad = sigi2(probe_TIP3P%atom)
      probe_TIP3P%q = q_seaton(probe_TIP3P%atom)
! Si(OH)4
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO,bondl_OH,angle_SiOH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigi2(probe_SiOH4%atom)
      probe_SiOH4%q = q_seaton(probe_SiOH4%atom)
      probe_SiOH4%q(1) = qsi(4)
! SiO4, just the 4 OH Oxygens
      call new_probe_mol(probe_SiO4,4,sigma_OH*0.5_wp)
      call tetra_coords(bondl_SiO,probe_SiO4%r)
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigi2(probe_SiO4%atom)
! O2
      call new_probe_mol(probe_O2,2,sig_O_O2*0.5_wp)
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_O2*0.5_wp /)
! N2
      call new_probe_mol(probe_N2,2,sig_N_N2*0.5_wp)
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_N2*0.5_wp /)
! CO2
      call new_probe_mol(probe_CO2,3,sig_O_CO2*0.5_wp)
      probe_CO2%rad(1) = sig_C_CO2*0.5_wp
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp, -bondl_CO2 /)
! Si(OH)4
      write(*,*) 'SiOH4'
      write(*,*) probe_SiOH4%n
      write(*,'(9i6)') probe_SiOH4%atom
      write(*,'(9f12.6)') probe_SiOH4%rad
      write(*,'(9f12.6)') probe_SiOH4%q
      write(*,'(3f12.6)') probe_SiOH4%r
! H2O
      write(*,*) 'H2O TIP3P'
      write(*,*) probe_TIP3P%n
      write(*,'(9i6)') probe_TIP3P%atom
      write(*,'(9f12.6)') probe_TIP3P%rad
      write(*,'(9f12.6)') probe_TIP3P%q
      write(*,'(3f12.6)') probe_TIP3P%r
! SiO4
      write(*,*) 'SiO4'
      write(*,*) probe_SiO4%n
      write(*,'(9f12.6)') probe_SiO4%rad
      write(*,'(3f12.6)') probe_SiO4%r
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

   SUBROUTINE WRITE_PROBE_MOL(iu,P)
      USE atom_types
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      integer:: i
      write(iu,*) P%n
      write(iu,*) 'Molecule'
      do i = 1,P%n
!        write(iu,'(a6,i5,a4,2x,a3,i6,4x,3f8.3')'HETATM',i,atom_name(P%atom(i)), &
!                                               '   ',i,P%r(1:3,i)*angstrom
         write(iu,'(a2,3(1x,f14.8))') atom_name(P%atom(i)),P%r(1:3,i)*Angstrom
      end do
   END SUBROUTINE

END MODULE probe_mol_mod

