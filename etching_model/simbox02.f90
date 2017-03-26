
include 'precision_mod.f90'
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
      integer,parameter:: natom_unit=24
      integer,parameter:: nunitc = 4
      real(wp):: r3(3)
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:),proximity2(:,:)
      real(wp):: CL,bl1,bl2,r3l
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      type(atom_list):: lst_positive_x,lst_negative_x,lstSi
      type(atom_list):: lst_positive_y,lst_negative_y
      type(atom_list):: lst_positive_z,lst_negative_z
!
      open(unit=16,file='connect_struct.out')
      open(unit=14,file='simbox.out')
      open(unit=15,file='neighbor.out')
      n_silicon_etch = 320
      natom_max = natom_unit*nunitc*nunitc*nunitc
      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(proximity(natom_max,4),proximity2(natom_max,4))
      proximity = 0
      write (*,*) natom_max

! read input data
      open(unit=13,FILE = 'unitcell.dat')
      do i = 1,natom_unit
         read(13,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         if (i <= 8) then
            atom(i) = iSilicon
         else
            atom(i) = iOxygen
         end if
      end do
      close(13,status='keep')

!------------------------------------------------------------
      CL = 40.0_wp
      boxl = nunitc*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      bl1 = CL*sqrt(2.0_wp)/8.0_wp
      bl2 = CL*sqrt(1.5_wp)/4.0_wp
!
! Creating the simulation box
!
      natom = 0
      do k = 1,nunitc
      do j = 1,nunitc
      do i = 1,nunitc
      do m = 1,natom_unit
         natom = natom + 1
         atom(natom) = atom(m)
         rxyz(natom,1) = rxyz(m,1) + (i - 1)*CL
         rxyz(natom,2) = rxyz(m,2) + (j - 1)*CL
         rxyz(natom,3) = rxyz(m,3) + (k - 1)*CL
      end do
      end do
      end do
      end do
      write (*,*) natom
! origin at centre of box
      rxyz(1:natom,1:3) = rxyz(1:natom,1:3) - CL*0.5_wp*(nunitc - 1)
      forall(i = 1:natom) rxyz(i,1:3) = rxyz(i,1:3) - boxl*anint(rxyz(i,1:3)*boxli)
      call write_xyz(14,0)
      close(14,status='KEEP')

! determine connectivity
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
                n = n + 1
                proximity(k,n) = j
            end if
         end do
      end do
      do i = 1,natom
         write (15,*) atom_name(atom(i)),proximity(i,1:4)
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

! Make array of silicon atoms
      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == iSilicon) call add_to_list(lstSi,i)
      end do
      print *,'no. of silicons = ',lstSi%n
      print *,count(atom == iSilicon)
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

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if


! determine connectivity ignoring pbc
      proximity2 = 0
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
                n = n + 1
                proximity2(k,n) = j
            end if
         end do
      end do

!
! Make list of bonds
      call bond_list

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
      if ((n_crossbond_x == 0).or. &
         (n_crossbond_y == 0).or. &
         (n_crossbond_z == 0)) then
         print *,'try again'
         stop
      end if

!
      do i = 1,n_crossbond_x
         ib = crossbond_x(i)
         if (rxyz(ibond(1,ib),1) > 0.0) then
            call add_to_list(lst_positive_x,ibond(1,ib))
            call add_to_list(lst_negative_x,ibond(2,ib))
            if (rxyz(ibond(2,ib),1) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),1),rxyz(ibond(1,ib),1)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_x,ibond(1,ib))
            call add_to_list(lst_positive_x,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_y
         ib = crossbond_y(i)
         if (rxyz(ibond(1,ib),2) > 0.0) then
            call add_to_list(lst_positive_y,ibond(1,ib))
            call add_to_list(lst_negative_y,ibond(2,ib))
            if (rxyz(ibond(2,ib),2) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),2),rxyz(ibond(1,ib),2)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_y,ibond(1,ib))
            call add_to_list(lst_positive_y,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_z
         ib = crossbond_z(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_positive_z,ibond(1,ib))
            call add_to_list(lst_negative_z,ibond(2,ib))
            if (rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_z,ibond(1,ib))
            call add_to_list(lst_positive_z,ibond(2,ib))
         end if
      end do
!

!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

!
! Finding a cluster which crosses all faces
!
      ic = 0
      do ib = 1,n_crossbond_x

         connected_x = .false.
         if (atomL(lst_positive_x%i(ib)) == atomL(lst_negative_x%i(ib))) then
            connected_x = .true.
         end if
         if (.not.connected_x) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(ib)) == atomL(lst_negative_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.not.connected_y) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(ib)) == atomL(lst_positive_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.not.connected_y) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(ib)) == atomL(lst_negative_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.not.connected_z) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(ib)) == atomL(lst_positive_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.not.connected_z) CYCLE

         spanning_cluster = connected_x.AND.connected_y.AND.connected_z
         if (spanning_cluster) then
            ic = atomL(lst_positive_x%i(ib))
            EXIT
         end if
      end do
      if (.not.spanning_cluster) then
         stop 'no spanning cluster found'
      end if




      if (spanning_cluster) then
         call new_frame_file(imve,'frame',1)
         call write_frame(imve,1,natom)
         natomc = 0
         do j = 1,natom
         if (atomL(j) == ic) then
            natomc = natomc + 1
            !write(16,*) atom_name(atom(j)),rxyz(j,1),rxyz(j,2),rxyz(j,3)
            !write(77,110)'HETATM',j,atom_name(atom(j)),'   ',j,rxyz(j,1:3)*angstrom
110         format(a6,i5,a4,2x,a3,i6,4x,3f8.3)
         end if
         end do
      end if

      if (.not.spanning_cluster) stop 'no spanning cluster found'

      if (spanning_cluster) then
         allocate(rxyzc(natomc,3),atomc(natomc))
         natomc = 0
         do j = 1,natom
            if (atomL(j) == ic) then
               natomc = natomc + 1
               rxyzc(natomc,1:3) = rxyz(j,:)
               atomc(natomc) = atom(j)
            end if
         end do
         atom(natomc + 1:) = 0
      end if
      deallocate(rxyz)
      deallocate(atom)
      allocate(rxyz(natomc*27,3),atom(natomc*27))
      rxyz(1:natomc,1:3) = rxyzc(1:natomc,1:3)
      atom(1:natomc) = atomc(1:natomc)


      natom = natomc
      proximity = 0
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)

            if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
               n = n + 1
               proximity(k,n) = j
            end if
         end do
      end do
      call new_frame_file(imve,'frame',2)
      call write_frame(imve,1,natom)
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster



      natom = natomc
      proximity = 0
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
               n = n + 1
               proximity(k,n) = j
            end if
         end do
      end do
      call new_frame_file(imve,'frame',2)
      call write_frame(imve,1,natom)




      ib = 0
      print *,'natomc = ',natomc
      do iz = -1,1
      do iy = -1,1
      do ix = -1,1
         if ((ix == 0).and.(iy == 0).and.(iz == 0)) cycle
         ib = ib + 1
print *,ib,ib*natomc + 1,(ib + 1)*natomc
         rxyz(ib*natomc + 1:(ib + 1)*natomc,1) = rxyzc(1:natomc,1) + boxl*ix
         rxyz(ib*natomc + 1:(ib + 1)*natomc,2) = rxyzc(1:natomc,2) + boxl*iy
         rxyz(ib*natomc + 1:(ib + 1)*natomc,3) = rxyzc(1:natomc,3) + boxl*iz
         atom(ib*natomc + 1:(ib + 1)*natomc) = atomc(1:natomc)
      end do
      end do
      end do
      natom = 27*natomc
      print *,'natom = ',natom

      open(99,file='connected_27.xyz')
      call write_xyz(99,0)

      boxl = boxl*3
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl

      deallocate(proximity)
      allocate(proximity(natom,4))
      proximity = 0
      do k = 1,natom  ! - 1
         n = 0
         do j = 1,natom      !k + 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            r3l = dot_product(r3,r3)
            if ( (r3l <= bl1*bl1).OR.(r3l <= bl2*bl2) ) then
                n = n + 1
                proximity(k,n) = j
            end if
         end do
      end do

      deallocate(nbond,ibond)
      allocate( ibond(2,natom*4),nbond(natom) )
      call new_frame_file(imve,'frame',3)
      call write_frame(imve,1,natom)

      CLOSE(UNIT = 16,STATUS = 'KEEP')

CONTAINS

   pure FUNCTION len_3d(r)
      real(wp):: len_3d
      real(wp),intent(in):: r(3)
        len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
   END FUNCTION

END PROGRAM SIMBOX

