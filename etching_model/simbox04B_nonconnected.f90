
include 'precision_mod.f90'
include 'command_line_mod.f90'
include 'command_arg_mod.f90'
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
      USE command_arg_mod
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
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: NAvo = 6.02214e23_wp
      real(wp):: n_si,n_ob,n_oh
      real(wp):: r3(3),rn1,rn2,rn3
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      integer:: natom_old
      integer,allocatable:: atom_old(:),proximity_old(:,:),proximity2(:,:)
      integer,allocatable:: atom_old1(:)
      real,allocatable:: rxyz_old(:,:)
      integer,allocatable:: atomc(:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l,frem,dens_gcm3,volbox,fdens_wt
      integer:: narg,len,status,noxdel
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz,isil
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl,nc
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz,irand_seed0
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      character(len=132) infile,ctmp,outfile,c6*6,c5*5,ftyp*3
      type(atom_list):: lstSi
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_com_arg(0, ctmp, "command")
      if (narg /= 4) then
         write(*,*)'usage :'
         write(*,*) trim(ctmp),' frem irand_seed infile outfile'
         write(*,*) 'frem = fractional density (weight) required'
         write(*,*) 'irand_seed = seed for RNG'
         write(*,*) 'infile = name of input file'
         write(*,*) 'outfile = base name of output file'
         stop 'error: incorrect no. of command line arguments'
      end if
      call get_com_arg(1, frem, "fractional density (weight) required")
      call get_com_arg(2, irand_seed0, "seed for RNG")
      call get_com_arg(3, infile, "input file")
      call get_com_arg(4, outfile, "base name of output file")
      irand_seed = irand_seed0
      write(*,'(a,2f9.6)') 'Check: first 2 random numbers are ',rand(),rand()
!
      open(unit=14,file=trim(infile))
      nseed = 0
!
      nc = len_trim(infile)
print *,"'",trim(infile),"'"
print *,nc
      ftyp = infile(nc - 2:nc)
      select case(ftyp)
      case('xyz')
         read (14,*) natom
         read (14,*) boxl
         open(unit=15,file=trim(infile(1:nc - 3))//'con')
      case('pdb')
         read (14,*) c6,natom  ;print *,'natom = ',natom
         read (14,*) c6,boxl   ;print *,' boxl = ',boxl
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
!
      natom_max = natom
      boxl = boxl/Angstrom
      boxl2 = boxl*0.5_wp
      boxli = 1.0_wp/boxl
      volbox = (boxl*boxl*boxl)*(Angstrom**3)
!
      allocate(rxyz(1:natom_max,3),rxyz_old(1:natom_max,3))
      allocate(atom(1:natom_max),atom_old(1:natom_max))
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(cluster_list_x(natom_max),cluster_list_y(natom_max),cluster_list_z(natom_max))
      allocate(proximity(natom_max,4),proximity_old(natom_max,4),proximity2(natom_max,4))
!
      select case(ftyp)
      case('xyz')
         do i = 1,natom
            read (14,*) ctmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
         close(14)
      case('pdb')
         do i = 1,natom
            read (14,'(a6)',advance = 'no') c6
            read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
      end select
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom

      proximity = 0
      select case(ftyp)
      case('xyz')
         read(15,*) itmp
         do i = 1,natom
            read(15,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(15)
      case('pdb')
         do i = 1,natom
            read(14,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(14)
      end select
!
! randomly translate the initial structure
!     rn1 = rand()
!     rn2 = rand()
!     rn3 = rand()
!     rxyz(1:natom,1) = rxyz(1:natom,1) + (2.0_wp*rn1 - 1.0_wp)*boxl
!     rxyz(1:natom,2) = rxyz(1:natom,2) + (2.0_wp*rn2 - 1.0_wp)*boxl
!     rxyz(1:natom,3) = rxyz(1:natom,3) + (2.0_wp*rn3 - 1.0_wp)*boxl
!     call pbc(rxyz)
!
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      write(*,*) 'density % = ',100.0*frac_dens_wt()
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
!
! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
!      if (n_cluster /= 1) then
!         stop 'initially there should be one cluster'
!      end if

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

!do i = 1,natom
!if (atom(i) ==iSilicon) then
!   do j =1,4
!      ii = proximity(i,j)
!      if (ii == 0) STOP '0 error in proximity'
!   end do
!end if
!end do
!do i = 1,natom
!if ( all(proximity(i,:) == 0) ) then
!   write(*,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
!   stop 'error in proximity'
!end if
!end do

! Make array of silicon atoms
      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == isilicon) call add_to_list(lstSi,i)
      end do
      print *,'no. of silicons = ',lstSi%n
      print *,count(atom(1:natom) == isilicon)
!
! Perfom the etching procedure
!
      isil = 0
      etching: do
         isil = isil + 1
         ! pick a Silicon to remove
!        call rand_from_list(lstSi,iat)
         do
            iat = int(rand()*natom) + 1
            if (atom(iat) == iSilicon) exit
         end do
! print *,'delete atom ',iat
! store old config
         call store_config()
         call delete_atom(iat)
! remove the free oxygens
         j = natom
         noxdel = 0
         do
            if (atom(j) == iOxygen .or. atom(j) == iOxygenH) then
            if (count(proximity(j,:) > 0) < 1) then
               call delete_atom(j)
               noxdel = noxdel + 1
            end if
            end if
            j = j - 1
            if (j == 0) exit
         end do
! convert singly bonded oxygens to OH
         do i = 1,natom
            if (atom(i) == iOxygen) then
            if (count(proximity(i,:) > 0) < 2) then
               atom(i) = iOxygenH
            end if
            end if
         end do

!! check the clusters
!         atomL = 0
!         call Init_HKNonLattice(natom)
!         call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!         !print *,'n_cluster = ',n_cluster
!         if (n_cluster > 1) then
!!           call remove_from_list(lstSi,iat)
!!           if (atom(iat) ==iSilicon) then
!!              where(lstSi%i == (natom+1)) lstSi%i = iat
!!           end if
!            call restore_config()
!            cycle etching
!         end if

         ! check for errors
         if (ANY(proximity(1:natom,:) > natom)) then
            print *,'error: proximity(1:natom,:) > natom '
            do i = 1,natom
               if (ANY(proximity(i,:) > natom)) then
               print *,i,' : ',proximity(i,:)
               end if
            end do
            stop
         end if
         if (ANY(proximity(natom + 1:,:) /= 0)) then
            print *,'error: proximity(natom + 1:,:) /= 0 '
            stop
         end if

         fdens_wt = frac_dens_wt()
         if (mod(isil,10) == 0) print '(i6,f9.5)',isil,fdens_wt

         if (fdens_wt <= frem) exit etching

      end do etching
!
!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster


! check everything is ok
! check each Si has 4 bonds
      do i = 1,natom
         if ((atom(i) == iSilicon)) then
         do j = 1,4
            ii = proximity(i,j)
            if (ii == 0) STOP '5 error in proximity'
         end do
         end if
      end do
      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
      call check_proximity2(ok,ii)
      if (.not.ok) then
         call write_proximity(6,1)
         print *,'proximity array is inconsistent (2) ',ii
         stop
      end if
!
! calculate the density
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      write(*,*) 'density % = ',100.0*frac_dens_wt()
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
      write(unit=c5,fmt='(f5.3)') frem
      open(unit=imve,file=trim(outfile)//'_'//c5//'_'//trim(c6)//'_'//trim(ctmp)//'.pdb')
      call write_frame(imve,1,natom)
      STOP 'normal termination'

CONTAINS

   FUNCTION frac_dens_wt()
      real(wp):: frac_dens_wt
      real(wp):: mass
      integer:: i
      mass = 0.0_wp
      do i = 1,natom
         select case(atom(i))
         case(iSilicon)
            mass = mass + mSi
         case(iOxygen)
            mass = mass + mOx
         case(iOxygenH)
            mass = mass + mOH
         case default
            stop 'unknown atom type'
         end select
      end do
      frac_dens_wt = mass/(8000*mSi + 16000*mOx)
   END FUNCTION frac_dens_wt

   SUBROUTINE store_config()
      atom_old = atom
      rxyz_old = rxyz
      proximity_old = proximity
      natom_old = natom
   END SUBROUTINE

   SUBROUTINE restore_config()
      atom = atom_old
      rxyz = rxyz_old
      proximity = proximity_old
      natom = natom_old
   END SUBROUTINE

END PROGRAM SIMBOX

