
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


MODULE condensation_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE condensation(iO1,iO2)
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      integer,intent(inout):: iO1,iO2
      integer:: i,ii,ilst(3),itmp
      integer:: iSi2,iSi1
!
!      H          H        O
!      |    +     |  -->  / \  + H2O
!   Si-O       Si-O
!
!
!      |        |           iO1           iO2
!     iO1   +  iO2   -->    /  \   +      /  \
!     /        /           /    \
!   iSi1     iSi2        iSi1  iSi2
!
!! iO2 cannot be a seed
!      if (iO2 <= nseed) then
!         if (iO1 <= nseed) RETURN  ! don't allow 2 seed OH groups to react
!         itmp = iO2
!         iO2 = iO1
!         iO1 = itmp
!      end if

!     Pick iSi1, the Silicon attached to iO1
      iSi1 = 0
      do i = 1,ncmax(atom(iO1))
         ii = proximity(iO1,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi1 = ii
         end if
      end do
      iSi2 = 0
      do i = 1,ncmax(atom(iO2))
         ii = proximity(iO2,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi2 = ii
         end if
      end do
!
      atom(iO1) = iOxygen
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO1, 0, iSi2)
      call delete_atom(iO2)
      RETURN
   END SUBROUTINE

END MODULE condensation_mod


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
      USE condensation_mod
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
      real,allocatable:: rxyz_old(:,:)
      integer,allocatable:: atomc(:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l,frem,dens_gcm3,volbox,fdens_wt
      real(wp):: r1(3),r2(3),rmax,rmax2,rl1,rl2,rmaxSiSi,rmaxSiSi2
      integer:: narg,len,status,noxdel,iSi1,iSi2
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz,isil,jj
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl,lst4(4),nc
      integer:: ipairs(4),ip(3,4),itry(3),ntry,ia,i1,i2,a1,a2
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz,irand_seed0
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster,success
      character(len=132) infile,ctmp,outfile,c6*6,c5*5,ftyp*3
      type(atom_list):: lstSi
!
      ip(1,:) = (/ 1,2,3,4 /)
      ip(2,:) = (/ 1,4,3,2 /)
      ip(3,:) = (/ 1,3,2,4 /)
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_com_arg(0, ctmp, "command")
      if (narg /= 6) then
         write(*,*)'usage :'
         write(*,*) trim(ctmp),' frem rmax rmaxSiSi irand_seed  infile  outfile'
         write(*,*)
         write(*,*) 'frem = fractional density (weight) required'
         write(*,*) "rmax = maximum radius for condensing OH's [Angstrom]"
         write(*,*) 'max radius for Silicas [Angstrom]'
         write(*,*) 'irand_seed = seed for RNG'
         write(*,*) 'infile = name of input file'
         write(*,*) 'outfile = base name of output file'
         stop 'error: incorrect no. of command line arguments'
      end if
      call get_com_arg(1, frem, "fractional density (weight) required")
      call get_com_arg(2, rmax, "max radius for condensing OH's")
      call get_com_arg(3, rmaxSiSi, "max radius for Silicas")
      call get_com_arg(4, irand_seed0, "seed for RNG")
      call get_com_arg(5, infile, "input file")
      call get_com_arg(6, outfile, "base name of output file")
      irand_seed = irand_seed0
      write(*,'(a,2f9.6)') 'Check: first 2 random numbers are ',rand(),rand()
      rmax = rmax/Angstrom
      rmax2 = rmax**2
      rmaxSiSi = rmaxSiSi/Angstrom
      rmaxSiSi2 = rmaxSiSi**2

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

!open(29,file='neighbor.out')
!do i = 1,natom
!   read(29,*) ctmp,(proximity(i,j),j = 1,4)
!   atom(i) = name2atom(trim(ctmp))
!end do
!do i = 1,natom
!   read(14,*) ctmp,rxyz(i,:)
!end do
!  rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
!
!      rn1 = rand()
!      rn2 = rand()
!      rn3 = rand()
!      rxyz(1:natom,1) = rxyz(1:natom,1) + (2.0_wp*rn1 - 1.0_wp)*boxl
!      rxyz(1:natom,2) = rxyz(1:natom,2) + (2.0_wp*rn2 - 1.0_wp)*boxl
!      rxyz(1:natom,3) = rxyz(1:natom,3) + (2.0_wp*rn3 - 1.0_wp)*boxl
!      call pbc(rxyz)
!
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(1000*mSi + 2000*mOx)
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

         ! Make array of silicon atoms
         lstSi%n = 0
         do i = 1,natom
            if (atom(i) == isilicon) call add_to_list(lstSi,i)
         end do

         ! pick a Silicon to remove
100      CONTINUE
         ! print *,'lstSi%n ',lstSi%n
         if (lstSi%n == 0) STOP 'NO suitable Silicon found for removal'
         call rand_from_list(lstSi,iat)
!do
!   iat = int(rand()*natom) + 1
!   if (atom(iat) == iSilicon) exit
!end do
         lst4 = proximity(iat,:)

! store old config
         call store_config()

         where (lst4 == natom) lst4 = iat
         call delete_atom(iat)

         itry = (/ 1,2,3 /)
         ntry = 3
         success = .false.
         do ii = 1,3
            jj = int(rand()*ntry) + 1
            ipairs = lst4(ip(itry(jj),:))
            RING_if: if (OH_groups_OK(ipairs(1),ipairs(2)).and.OH_groups_OK(ipairs(3),ipairs(4))) then
               ! Check if the 2 pairs of OH are are close enough
               r1 = rxyz(ipairs(1),:) - rxyz(ipairs(2),:)
               r2 = rxyz(ipairs(3),:) - rxyz(ipairs(4),:)
               call pbc(r1)
               call pbc(r2)
               rl1 = dot_product(r1,r1)
               rl2 = dot_product(r2,r2)
               OH_if: if ( max(rl1,rl2) < rmax2 ) then
                  ! Check if the 1st two Si are are close enough
                  iSi1 = 0
                  do i = 1,2
                     k = proximity(ipairs(1),i)
                     if (k == 0) cycle
                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                     iSi1 = k
                  end do
                  iSi2 = 0
                  do i = 1,2
                     k = proximity(ipairs(2),i)
                     if (k == 0) cycle
                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                     iSi2 = k
                  end do
                  r1 = rxyz(iSi1,:) - rxyz(iSi2,:)
                  call pbc(r1)
                  rl1 = dot_product(r1,r1)
                  Si_1_if: if ( rl1 < rmaxSiSi2 ) then
                     ! Check if the 2nd two Si are are close enough
                     iSi1 = 0
                     do i = 1,2
                        k = proximity(ipairs(3),i)
                        if (k == 0) cycle
                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                        iSi1 = k
                     end do
                     iSi2 = 0
                     do i = 1,2
                        k = proximity(ipairs(4),i)
                        if (k == 0) cycle
                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                        iSi2 = k
                     end do
                     r2 = rxyz(iSi1,:) - rxyz(iSi2,:)
                     call pbc(r2)
                     rl2 = dot_product(r2,r2)
                     Si_2_if: if ( rl2 < rmaxSiSi2 ) then
                        success = .true.
                        exit
                     else Si_2_if
                        ! print *,'Si - Si 2 failed'
                     end if Si_2_if
                  else Si_1_if
                     ! print *,'Si - Si 1 failed'
                  end if Si_1_if
               else OH_if
                  ! print *,'OH constraint failed'
               end if OH_if
            else RING_if
               ! print *,'Ring constraint failed'
               itry(jj) = itry(ntry)
               ntry = ntry - 1
            end if RING_if
         end do
         if (.not.success) then
            call restore_config()
            call remove_from_list(lstSi,iat)
            GOTO 100
!            print *,'could not find 2 suitable pairs to condense'
!            cycle etching
         end if

         r1 = rxyz(ipairs(2),:) - rxyz(ipairs(1),:)
         call pbc(r1)
         rxyz(ipairs(1),:) = rxyz(ipairs(1),:) +  r1*0.5_wp
         call pbc(rxyz(ipairs(1),:))

         r2 = rxyz(ipairs(4),:) - rxyz(ipairs(3),:)
         call pbc(r2)
         rxyz(ipairs(3),:) = rxyz(ipairs(3),:) +  r2*0.5_wp
         call pbc(rxyz(ipairs(3),:))

         if (maxval(ipairs(1:2)) > maxval(ipairs(3:4))) then
            call condensation(ipairs(1),ipairs(2))
            call condensation(ipairs(3),ipairs(4))
         else
            call condensation(ipairs(3),ipairs(4))
            call condensation(ipairs(1),ipairs(2))
         end if


!call check_proximity(ok,ii)
!if (.not.ok) then
!print *,'just after condensation 2'
!print *,'proximity array is inconsistent ',ii
!stop
!end if
!do i = 1,natom
!ia = atom(i)
!select case(ia)
!case(iSilicon)  ! Check Si is bonded to 4 atoms
!   do j = 1,ncmax(ia)
!      k = proximity(i,j)
!      if (k == 0) STOP 'ERROR 01'
!      if (atom(k) ==iSilicon) STOP 'ERROR 02'
!      if (atom(k) ==iHydrogen) STOP 'ERROR 03'
!   end do
!case(iOxygen)  ! Bridging Oxygen must be bonded to Si
!   i1 = proximity(i,1)
!   if (i1 == 0) then
!      print *,'i = ',i
!      print *,'ipairs = ',ipairs
!      STOP 'ERROR 04 B'
!   end if
!   a1 = atom(proximity(i,1))
!   if (a1 /= iSilicon) STOP 'ERROR 05'
!      i2 = proximity(i,2)  !  must be bonded to  2 Silicons
!      if (i2 == 0) STOP 'ERROR 06'
!      a2 = atom(proximity(i,2))
!      if (a2 /= iSilicon) STOP 'ERROR 07'
!   if (any(proximity(i,3:) /= 0)) STOP 'ERROR 08'
!case default
!   stop 'unknown atom type'
!end select
!end do


!! remove the free oxygens
!         j = natom
!         noxdel = 0
!         do
!            if (atom(j) == iOxygen .or. atom(j) == iOxygenH) then
!            if (count(proximity(j,:) > 0) < 1) then
!               call delete_atom(j)
!               noxdel = noxdel + 1
!            end if
!            end if
!            j = j - 1
!            if (j == 0) exit
!         end do
!! convert singly bonded oxygens to OH
!         do i = 1,natom
!            if (atom(i) == iOxygen) then
!            if (count(proximity(i,:) > 0) < 2) then
!               atom(i) = iOxygenH
!            end if
!            end if
!         end do

! check the clusters
         atomL = 0
         call Init_HKNonLattice(natom)
         call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
         !print *,'n_cluster = ',n_cluster
         if (n_cluster > 1) then
!           call remove_from_list(lstSi,iat)
!           if (atom(iat) ==iSilicon) then
!              where(lstSi%i == (natom+1)) lstSi%i = iat
!           end if
            call restore_config()
            cycle etching
         end if

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

      call bond_list
      proximity2 = proximity
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (dot_product(r3,r3) > boxl2*boxl2 ) then
            do j = 1,4
               if (proximity2(ibond(1,i),j) == ibond(2,i)) then
                   proximity2(ibond(1,i),j) = 0
               end if
            end do
            do j = 1,4
               if (proximity2(ibond(2,i),j) == ibond(1,i)) then
                   proximity2(ibond(2,i),j) = 0
               end if
            end do
         end if
         if (ABS(r3(1)) > boxl2) then
            n_crossbond_x = n_crossbond_x + 1
            crossbond_x(n_crossbond_x) = i
         end if
         if (ABS(r3(2)) > boxl2) then
            n_crossbond_y = n_crossbond_y + 1
            crossbond_y(n_crossbond_y) = i
         end if
         if (ABS(r3(3)) > boxl2) then
            n_crossbond_z = n_crossbond_z + 1
            crossbond_z(n_crossbond_z) = i
         end if
      end do
      atomL = 0

      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity2(1:natom,:),n_cluster,atomL)

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
!
! Finding a cluster which crosses all faces
!
       icx_n = 0
       icy_n = 0
       icz_n = 0
       do ibx = 1,n_crossbond_x
          icx = crossbond_x(ibx)
          if ( atomL(ibond(1,icx)) == atomL(ibond(2,icx)) ) then
             icx_n = icx_n + 1
             cluster_list_x(icx_n) = atomL(ibond(1,icx))
          end if
       end do

       do iby = 1,n_crossbond_y
          icy = crossbond_y(iby)
          if ( atomL(ibond(1,icy)) == atomL(ibond(2,icy)) ) then
             icy_n = icy_n + 1
             cluster_list_y(icy_n) = atomL(ibond(1,icy))
         end if
       end do

       do ibz = 1,n_crossbond_z
          icz = crossbond_z(ibz)
          if ( atomL(ibond(1,icz)) == atomL(ibond(2,icz)) ) then
             icz_n = icz_n + 1
             cluster_list_z(icz_n) = atomL(ibond(1,icz))
          end if
       end do

      spanning_cluster = .false.
      outer: do icl = 1, icx_n
      do icly = 1, icy_n
      do iclz = 1, icz_n
         if ((cluster_list_x(icl) == cluster_list_y(icly)) .AND. &
             (cluster_list_x(icl) == cluster_list_z(iclz))) then
             spanning_cluster = .true.
             exit outer
         end if
      end do
      end do
      end do outer

      do icly = 1,natom
         if ( atomL(icly) == cluster_list_x(icl) ) exit
      end do

      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      j = natom
      do
         if (atomL(j) /= atomL(icly)) then
print *,'deleting atoms not in 1st spanning cluster'
            call delete_atom(j)
            atomL(j) = atomL(icly)
         end if
         j = j - 1
         if (j == 0) exit
      end do

      do i = 1,natom
         if ((atom(i) == iSilicon).and.(atomL(i) == atomL(icly))) then
         do j = 1,4
            ii = proximity(i,j)
            if (ii == 0) STOP '5 error in proximity'
         end do
         end if
      end do

      write(*,*) 'natom ',natom
      write(*,*) 'irand_seed ',irand_seed
      if (spanning_cluster) then
         imve = 77
         write(unit=ctmp,fmt='(i6.6)') irand_seed0
         write(unit=c6,fmt='(i6.6)') natom
         write(unit=c5,fmt='(f5.3)') frem
         open(unit=imve,file=trim(outfile)//'_'//c5//'_'//trim(c6)//'_'//trim(ctmp)//'.pdb')
         call write_frame(imve,1,natom)
      else
         print *,'A spanning cluster was NOT generated'
      end if


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
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(1000*mSi + 2000*mOx)
      write(*,*) 'density % = ',100.0*frac_dens_wt()
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
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
      frac_dens_wt = mass/(1000*mSi + 2000*mOx)
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

