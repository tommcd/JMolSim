
include 'precision_mod.f90'
include 'global_vars.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'files_mod.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'charges2_mod.f90'
include 'matvec3d_mod.f90'
include 'rotate_axis_mod.f90'
include 'nlist_mod.f90'
include 'tetra_coords_mod.f90'
include 'quat2mat_mod.f90'
include 'ran_point_sphere_mod.f90'
include 'probe_mol_mod_DVF.f90'
include 'lj_el_mod.f90'
include 'Henrys_law_calc_mod.f90'
include 'voidage4_mod.f90'
include 'deposition4_mod_H.f90'


PROGRAM SIMBOX
      USE precision_mod
      USE global_vars
      USE constants_mod
      USE atom_types
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      USE files_mod
      USE probe_mol_mod
      USE seaton_mod
      USE rand_mod
      USE atom_types
      USE lj_el_mod
      USE nlist_mod
      USE charges_mod
      USE seaton_mod
      USE charges_mod
      USE lj_el_mod
      USE voidage_mod
      USE Henrys_law_calc_mod
      USE deposition_mod
      implicit none
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      integer:: j,i,ip,k,itmp,narg,len,status
      real(wp),allocatable:: coprxyz(:,:)
      real(wp):: Temp_K,dlayer,x,y,z
      character:: c6*6,c5*5
      character(len=132) infile,ctmp
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', ctmp(1:len)
      if(narg < 3)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' natom atom_file T[K]'
         stop
      end if
!      
      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) natom
      print *,'natom = ',natom
!      
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)      
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) Temp_K
!
      kboltzT = K_ev*Temp_K
      print *,'kT [eV] = ',kboltzT
!
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      natom_max = natom
      write (*,*) natom_max
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(atom(1:natom_max))
      allocate(charge(natom))
      allocate(proximity(natom_max,4))

! Read in coordinates and connectivity from pdb file
      do i=1,natom
         read (14,'(a6)',advance='no') c6
         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j=1,3)
         atom(i) = name2atom(trim(ctmp))             
      end do
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
      proximity = 0
      do i=1,natom
         read(14,'(a32)')ctmp
         do j = 1,ncmax(atom(i))
            c5 = ctmp(6+5*(j)+1:6+5*(j)+5)
            read( unit=c5,fmt=* ) proximity(i,j)
         end do
      end do

      call assign_charge(1,natom)
      
write(*,*) 'sum charge = ',sum(charge(1:natom))
write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
write(*,*) 'n Si = ',count(atom==iSilicon)
write(*,*) 'n O  = ',count(atom==iOxygen)
write(*,*) 'n OH = ',count(atom==iOxygenH)
write(*,*) 'ntot = ',count(atom==iSilicon) + count(atom==iOxygen) + count(atom==iOxygenH)
!    
      call init_probe_mols()
      call INIT_NLIST(boxl,boxl,boxl,5.0_wp/angstrom)
      call NEW_NLIST
      call LJ_INIT
      dlayer = (5.0_wp/angstrom)
      call henry_profile(probe_SiO4,dlayer,1)
      call henry_profile(probe_tip3p,dlayer,2)


      open(unit=20,file='khenry.out')
      do ip = 1,2
      do k = 1,ncelx
         z = (k-0.5_wp)*delz
         do j = 1,ncely
            y = (j-0.5_wp)*dely
            do i = 1,ncelz
               x = (i-0.5_wp)*delx
               write(20,*) x,y,z,henryk(i,j,k,ip)
            end do
         end do
      end do
      write(20,'(/)')
      end do

      call den_profile_rxn(1,natom,dlayer)

      close(14,status='KEEP')
contains

   function command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   end function
   
   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
      CHARACTER(LEN=1000) :: TMPVAL
      INTEGER :: IARGC
      EXTERNAL   IARGC
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      ELSE IF (NUMBER > IARGC()) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 2
          RETURN
      END IF
      IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
      IF (PRESENT(LENGTH)) THEN
          IF (PRESENT(VALUE)) THEN
              LENGTH = LEN_TRIM(VALUE)
          ELSE
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          END IF
      END IF
      IF (PRESENT(STATUS)) STATUS = 0
      RETURN
   END SUBROUTINE GET_COMMAND_ARGUMENT
end program simbox


