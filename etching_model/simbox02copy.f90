
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
      real(wp):: r3(3)
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      character(6):: c6,ctmp*32
      character(5):: c5
      integer,allocatable:: atomc(:),proximity2(:,:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      open(unit=14,file='simbox.out')
      open(unit=15,file='vink_8_3000.pdb')
      n_silicon_etch = 1660
write(*,'(a)',advance = 'no')'enter random seed:'
read(*,*) irand_seed
write(*,'(a,f8.6)') 'Check: first random number is ',rand()
!      natom_max = natom_unit*nunitc*nunitc*nunitc
       natom_max = 24000
      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(cluster_list_x(natom_max),cluster_list_y(natom_max),cluster_list_z(natom_max))
      allocate(proximity(natom_max,4),proximity2(natom_max,4))
      proximity = 0
      write (*,*) natom_max
!      CL = 0.7168_wp*0.162_wp/0.1552_wp
      boxl = 7.13286_wp
!      boxl = 8.0_wp*CL
      boxl2 = boxl/2
      boxli = 1.0_wp/boxl
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
! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
!         stop 'there should be one cluster'
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
!      print *,count(atom==iSilicon)
!
! Perfom the etching procedure
!
      do i = 1, n_silicon_etch
         ! pick a Silicon to remove
         !call rand_from_list(lstSi,iat)
         do
            iat = int(rand()*natom) + 1
            if (atom(iat) == iSilicon) exit
         end do
!print *,'delete atom ',iat
         call delete_atom(iat)

         !call remove_from_list(lstSi,iat)
         !if (atom(iat) == iSilicon) then
         !   where(lstSi%i == (natom + 1)) lstSi%i = iat
         !end if
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


!      do i = 1,natom
!         if ( all(proximity(i,:) == 0) ) then
!            write(*,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
!            stop '1error in proximity'
!         end if
!      end do
!


      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if

! Make list of bonds
      call bond_list
      proximity2 = proximity
!
! Make arrays of crossing bonds
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (len_3d(r3) > boxl2 ) then
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
!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity2(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster






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





! do i = 1,natom
!!    if ( atom(i) /= iSilicon ) cycle
!    if ( atomL(i) /= cluster_list_x(icl) ) cycle
!    do j =1,4
!!       ii = proximity(i,j)
!!       if (ii == 0) STOP '3A: error in proximity'
!       ii = proximity2(i,j)
!       if (ii == 0) cycle
!       if ( (atomL(proximity2(i,j)) /= cluster_list_x(icl))) STOP '3B: error in proximity'
!    end do
! end do
!stop '3B:no errors in proximity'

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


      do i = 1,natom
          if (atom(i) == iOxygen) then
              if (count(proximity(i,:) > 0) < 2) then
                 atom(i) = iOxygenH
              end if
          end if
      end do





      write(*,*) natom
      if (spanning_cluster) then
         call new_frame_file(imve,'frame',1)
         call write_frame(14,1,natom)
      end if

CONTAINS

   pure FUNCTION len_3d(r)
      real(wp):: len_3d
      real(wp),intent(in):: r(3)
        len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
   END FUNCTION

END PROGRAM SIMBOX

