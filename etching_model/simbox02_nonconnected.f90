
include 'precision_mod.f90'
include 'Keating_parameters_vonAlfthan.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'HKNonLattice2.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'bond_list_mod.f90'
include 'frames_mod.f90'


PROGRAM SIMBOX
      USE precision_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE HKNonLattice_mod
      USE frames_mod
      implicit none
!      integer,parameter:: natom_unit=24
!      integer,parameter:: nunitc = 8
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: NAvo = 6.02214e23_wp
      real(wp):: n_si,n_ob,n_oh
      real(wp):: r3(3)
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l,frem,dens_gcm3,volbox
      integer:: narg,len,status,noxdel
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz,irand_seed0
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      character(len=132) infile,ctmp,outfile,c6*6,c5*5
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', ctmp(1:len)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' frem  irand_seed  out_file'
         stop
      end if
!
      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) frem
      print *,'frem = ',frem
!
      call get_command_argument (2, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) irand_seed
      irand_seed0 = irand_seed
      print *,'irand_seed = ',irand_seed
      write(*,'(a,2f9.6)') 'Check: first 2 random numbers are ',rand(),rand()
!
      call get_command_argument (3, outfile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', outfile(1:len)
!
      open(unit=15,file='vink_8_3000.pdb')
!     natom_max = natom_unit*nunitc*nunitc*nunitc
      natom_max = 24000

      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(cluster_list_x(natom_max),cluster_list_y(natom_max),cluster_list_z(natom_max))
      allocate(proximity(natom_max,4))
      proximity = 0
      write (*,*) natom_max
!      CL = 0.7168_wp*0.162_wp/0.1552_wp
      boxl = 7.13286_wp
!      boxl = 8.0_wp*CL
      boxl2 = boxl/2
      boxli = 1.0_wp/boxl
      volbox = (boxl*boxl*boxl)*1000_wp
      natom = natom_max
!      call read_xyz(14)
     read(15,*) c6   ;print *,c6
      do i = 1,natom
         read (15,'(a6)',advance = 'no') c6
         read (15,*) itmp,ctmp,itmp,rxyz(i,:)
!         print'(i6,a,i6,3f12.6)',i,trim(ctmp),itmp,rxyz(i,:)
         atom(i) = name2atom(trim(ctmp))
      end do
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
      do i = 1,natom
         read(15,'(a32)') ctmp
!print *,trim(ctmp)
!print *,ncmax(atom(i))
         do j = 1,ncmax(atom(i))
            c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
            read( unit=c5,fmt=* ) proximity(i,j)
!            read( unit=ctmp(6+5*(j-1)+1:6+5*(j-1)+5) ) proximity(i,j)
         end do
!         write(777,'(a6,5i6)')'CONECT',(proximity(i,j),j=1,ncmax(atom(i)))
      end do

!open(29,file='neighbor.out')
!do i = 1,natom
!   read(29,*) ctmp,(proximity(i,j),j = 1,4)
!   atom(i) = name2atom(trim(ctmp))
!end do
!do i = 1,natom
!   read(14,*) ctmp,rxyz(i,:)
!end do
!  rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom

      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
      n_silicon_etch = (1.0_wp - frem)*count(atom(1:natom) == iSilicon)  !1660
      write(*,*) 'no. of silicons to remove = ',n_silicon_etch
!
! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
         stop 'initially there should be one cluster'
      end if

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
      call check_proximity2(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent (2) ',ii
         stop
      end if


do i = 1,natom
if (atom(i) == iSilicon) then
   do j = 1,4
      ii = proximity(i,j)
      if (ii == 0) STOP '0 error in proximity'
   end do
end if
end do
!do i = 1,natom
!if ( all(proximity(i,:) == 0) ) then
!   write(*,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
!   stop 'error in proximity'
!end if
!end do


! Make array of silicon atoms
!      lstSi%n = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) call add_to_list(lstSi,i)
!      end do
!      print *,'no. of silicons = ',lstSi%n
!      print *,count(atom(1:natom) ==iSilicon)
!
! Perfom the etching procedure
!
      do i = 1, n_silicon_etch
         ! pick a Silicon to remove
         do
            iat = int(rand()*natom) + 1
            if (atom(iat) == iSilicon) exit
         end do
!print *,'delete atom ',iat
         call delete_atom(iat)

         ! check for errors
         if (ANY(proximity(1:natom,:) > natom)) then
            print *,'error: proximity(1:natom,:) > natom ',i
            stop
         end if
         if (ANY(proximity(natom + 1:,:) /= 0)) then
            print *,'error: proximity(natom + 1:,:) /= 0 ',i
            stop
         end if
      end do
!
      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if

!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      j = natom
      noxdel = 0
      do
         if (atom(j) == iOxygen) then
         if (count(proximity(j,:) > 0) < 1) then
            call delete_atom(j)
            noxdel = noxdel + 1
         end if
         end if
         j = j - 1
         if (j == 0) exit
      end do
      print *, 'no. of oxygens deleted = ',noxdel

do i = 1,natom
   if ((atom(i) == iSilicon)) then
   do j = 1,4
      ii = proximity(i,j)
      if (ii == 0) STOP '5 error in proximity'
   end do
   end if
end do

      do i = 1,natom
          if (atom(i) == iOxygen) then
              if (count(proximity(i,:) > 0) < 2) then
                 atom(i) = iOxygenH
              end if
          end if
      end do

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
!      call check_proximity2(ok,ii)
!      if (.not.ok) then
!         call write_proximity(6,1)
!         print *,'proximity array is inconsistent (2) ',ii
!         stop
!      end if
      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster


      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'

      write(*,*) natom
      imve = 77
      write(unit=ctmp,fmt='(i6.6)') irand_seed0
      write(unit=c6,fmt='(i6.6)') natom
      open(unit=imve,file=trim(outfile)//'_'//trim(c6)//'_'//trim(ctmp)//'.pdb')
      call write_frame(imve,1,natom)
      STOP 'normal termination'

CONTAINS

   FUNCTION command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   END FUNCTION

   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(len=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
      CHARACTER(len=1000) :: TMPVAL
      INTEGER :: IARGC
      EXTERNAL   IARGC
      if (NUMBER < 0) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(STATUS)) STATUS = 1
          RETURN
      else if (NUMBER > IARGC()) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(STATUS)) STATUS = 2
          RETURN
      end if
      if (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
      if (PRESENT(LENGTH)) then
          if (PRESENT(VALUE)) then
              LENGTH = LEN_TRIM(VALUE)
          else
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          end if
      end if
      if (PRESENT(STATUS)) STATUS = 0
      RETURN
   END SUBROUTINE GET_COMMAND_ARGUMENT

END PROGRAM SIMBOX

