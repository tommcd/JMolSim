
include 'precision_mod.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'Keating_parameters_vonAlfthan.f90'
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
      integer,parameter:: natom_unit=9
      integer,parameter:: nunitc = 2
      real(wp):: r3(3)
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp):: CL,bl1,bl2,r3l
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      type(atom_list):: lst_pos_x,lst_neg_x,lstSi
      type(atom_list):: lst_pos_y,lst_neg_y
      type(atom_list):: lst_pos_z,lst_neg_z
      integer:: idir = 1
!
      open(unit=16,file='connect_struct.out')
      open(unit=14,file='simbox.xyz')
      open(unit=15,file='neighbor.out')
      n_silicon_etch = 320
      natom_max = natom_unit*nunitc
      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(proximity(natom_max,4))
      proximity = 0
      write (*,*) 'natom_max = ',natom_max

! read input data
!      open(unit=13,FILE='unitcell.dat')
!      do i=1,natom_unit
!         read(13,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
!            atom(i) = iOxygen
!         end if
!      end do
!      close(13,status='keep')

!------------------------------------------------------------
      CL = 6.0_wp
      boxl = CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!
! Creating the simulation box
!
      natom = 0
      do k = 0,0
      do j = -2,2,2
      do i = -2,2,2
         natom = natom + 1
print *,i,j,k,natom
         atom(natom) = iOxygen
         rxyz(natom,1) = i
         rxyz(natom,2) = j
         rxyz(natom,3) = k
      end do
      end do
      end do
      atom(5) = iSilicon
      write (*,*) 'natom = ',natom

      call write_xyz(14,0)
      close(14,status='KEEP')

! determine connectivity
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            if ( (len_3d(r3) <= 1.1*2.0) ) then
                n = n + 1
                proximity(k,n) = j
            end if
         end do
      end do
      do i = 1,natom
         write (*,*) i,proximity(i,1:4)
      end do
      call new_frame_file(imve,'frame',0)
      call write_frame(imve,1,natom)


! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
         stop 'there should be one cluster'
      end if

!
! Make list of bonds
      call bond_list
      do i = 1,nbondtot
         print *,i,'   ',ibond(1:2,i)
      end do

!
! Make arrays of crossing bonds
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
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

!
! Make arrays of interface atoms
!
      write(*,*) n_crossbond_x
      write(*,*) n_crossbond_y
      write(*,*) n_crossbond_z
!
      do i = 1,n_crossbond_x
         ib = crossbond_x(i)
         if (rxyz(ibond(1,ib),1) > 0.0) then
            call add_to_list(lst_pos_x,ibond(1,ib))
            call add_to_list(lst_neg_x,ibond(2,ib))
            if (rxyz(ibond(2,ib),1) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),1),rxyz(ibond(1,ib),1)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_neg_x,ibond(1,ib))
            call add_to_list(lst_pos_x,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_y
         ib = crossbond_y(i)
         if (rxyz(ibond(1,ib),2) > 0.0) then
            call add_to_list(lst_pos_y,ibond(1,ib))
            call add_to_list(lst_neg_y,ibond(2,ib))
            if (rxyz(ibond(2,ib),2) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),2),rxyz(ibond(1,ib),2)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_neg_y,ibond(1,ib))
            call add_to_list(lst_pos_y,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_z
         ib = crossbond_z(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_pos_z,ibond(1,ib))
            call add_to_list(lst_neg_z,ibond(2,ib))
            if (rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_neg_z,ibond(1,ib))
            call add_to_list(lst_pos_z,ibond(2,ib))
         end if
      end do
!

!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      print *,'crossbond_x ',crossbond_x(1:n_crossbond_x)
      print *,'crossbond_y ',crossbond_y(1:n_crossbond_y)
      print *,'crossbond_z ',crossbond_z(1:n_crossbond_z)
      print *,'lst_pos_x ',lst_pos_x%i(1:lst_pos_x%n)
      print *,'lst_neg_x ',lst_neg_x%i(1:lst_neg_x%n)
      print *,'lst_pos_y ',lst_pos_y%i(1:lst_pos_y%n)
      print *,'lst_neg_y ',lst_neg_y%i(1:lst_neg_y%n)
      print *,'lst_pos_z ',lst_pos_z%i(1:lst_pos_z%n)
      print *,'lst_neg_z ',lst_neg_z%i(1:lst_neg_z%n)




      rxyz(natom + 1:2*natom,:) = rxyz(1:natom,:)
      rxyz(natom + 1:2*natom,idir) = rxyz(natom + 1:2*natom,idir) + boxl
      atom(natom + 1:2*natom) = atom(1:natom)
      proximity(natom + 1:2*natom,:) = proximity(1:natom,:) + natom
      do i = 1,n_crossbond_x
         call set_proximity(lst_pos_x%i(i),lst_neg_x%i(i),lst_neg_x%i(i) + natom)
         call set_proximity(lst_neg_x%i(i),lst_pos_x%i(i),lst_pos_x%i(i) + natom)
      end do

      boxl = boxl + boxl
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      natom = natom + natom




      call new_frame_file(imve,'frame',1)
      call write_frame(imve,1,natom)

!!
!      call Init_HKNonLattice(natom)
!      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!      print *,'n_cluster = ',n_cluster

CONTAINS

   pure FUNCTION len_3d(r)
      real(wp):: len_3d
      real(wp),intent(in):: r(:)
      len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
   END FUNCTION

END PROGRAM SIMBOX

