
include 'precision_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'HKNonLattice.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'bond_list_mod.f90'


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
      implicit none
      integer,parameter:: natomt = 1536
      real(wp) rxyz1(natomt,3),r3(3)
      integer crossbond_x(natomt),crossbond_y(natomt),crossbond_z(natomt)
      real(wp) CL
      integer j,i,m,n,ib,k,iat,ifirst
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      logical:: connected_x,connected_y,connected_z,consistent
      type(atom_list) :: lst_positive_x,lst_negative_x,lstSi
      type(atom_list) :: lst_positive_y,lst_negative_y
      type(atom_list) :: lst_positive_z,lst_negative_z
      character(2):: ctmp

      OPEN(unit=13,FILE = 'unitcell.dat')
      OPEN(unit=14,FILE = 'simbox.out')
      OPEN(unit=16,FILE = 'connect_struct.out')
      OPEN(unit=15,FILE = 'neighbor.out')
      natom = natomt
      natom_max = natomt
      allocate(rxyz(1:natom,3))
      allocate(atom(1:natom))
      allocate(proximity(natom,4))
      proximity = 0

! READ 'input.dat'
      do i = 1,24
         READ(13,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         if (i <= 8) then
            atom(i) = iSilicon
         else
            atom(i) = iOxygen
         end if
      end do
      CLOSE(UNIT = 13,STATUS = 'KEEP')

!------------------------------------------------------------
      CL = 40.0d0
      boxl = 4.0d0*CL
      boxli = 1.0d0/boxl

!    *******************************
!***** Creating the simulation box *
!    *******************************
      do k = 1,4
      do j = 1,4
      do i = 1,4
      do m = 1,24
         write(14,*) atom_name(atom(m)), &
                     rxyz(m,1) + (i - 1)*CL, &
                     rxyz(m,2) + (j - 1)*CL, &
                     rxyz(m,3) + (k - 1)*CL
      end do
      end do
      end do
      end do
      CLOSE(UNIT = 14,STATUS = 'KEEP')

!    ********************************
!***** Making list of the neighbors *
!    ********************************

      OPEN(unit=14,FILE = 'simbox.out')
      do i = 1,natom
         READ(14,*) ctmp, rxyz(i,1), rxyz(i,2), rxyz(i,3)
         atom(i) = name2atom(ctmp)
         rxyz1(i,1) = rxyz(i,1) - 1.5d0*CL
         rxyz1(i,2) = rxyz(i,2) - 1.5d0*CL
         rxyz1(i,3) = rxyz(i,3) - 1.5d0*CL
!--------Periodic boundary conditions
         rxyz(i,1) = rxyz1(i,1) - boxl*anint(rxyz1(i,1)*boxli)
         rxyz(i,2) = rxyz1(i,2) - boxl*anint(rxyz1(i,2)*boxli)
         rxyz(i,3) = rxyz1(i,3) - boxl*anint(rxyz1(i,3)*boxli)
      end do
!
      write (*,*) natom

      do k = 1,natom
         n = 0
         do j = 1,natom
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            if (((((len_3d(r3) <= (CL*sqrt(2.0d0)/8.0d0)) &
               .OR.(len_3d(r3) <= (CL*sqrt(3.0d0/2.0d0)/4.0d0))) &
               .AND.(j /= k))).AND.(atom(k) /= atom(j))) then
                n = n + 1
                proximity(k,n) = j
            end if
         end do
      end do

      do i = 1,natom
         write (15,*) atom_name(atom(i)),proximity(i,1:4)
      end do

!    *********************************
!***** Making array of silicon atoms *
!    *********************************

      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == iSilicon) call add_to_list(lstSi,i)
      end do

!    ***************************************
!***** Perfomance of the etching procedure *
!    ***************************************

      n_silicon_etch = 300

      do i = 1, n_silicon_etch
         call rand_from_list(lstSi,iat)
         call delete_atom(iat)
         call remove_from_list(lstSi,iat)
      end do

call check_proximity(consistent,ifirst)
if (.NOT.consistent) then
print *,'ifirst = ',ifirst
stop
end if

!    ****************************
!***** Making list of the bonds *
!    ****************************
      call bond_list

!    ***********************************
!***** Making arrays of crossing bonds *
!    ***********************************
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (ABS(r3(1))>boxl/2.0d0) then
            n_crossbond_x = n_crossbond_x + 1
            crossbond_x(n_crossbond_x) = i
         end if
         if (ABS(r3(2))>boxl/2.0d0) then
            n_crossbond_y = n_crossbond_y + 1
            crossbond_y(n_crossbond_y) = i
         end if
         if (ABS(r3(3))>boxl/2.0d0) then
            n_crossbond_z = n_crossbond_z + 1
            crossbond_z(n_crossbond_z) = i
         end if
      end do

!    *************************************
!***** Making arraies of interface atoms *
!    *************************************
      write(*,*) n_crossbond_x
      write(*,*) n_crossbond_y
      write(*,*) n_crossbond_z
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

!     ******************************************
!****** Applying the Hoschen-Kopelman algoritm *
!     ******************************************

      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)


!     *************************************************
!****** Finding the cluster which contains all planes *
!     *************************************************

      do i = 1,n_crossbond_x

         connected_x = .false.
         do j = 1,n_crossbond_x
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_x%i(j))) then
               connected_x = .true.
               exit
            end if
         end do
         if (.NOT.connected_x) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.NOT.connected_y) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(i)) == atomL(lst_positive_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.NOT.connected_y) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.NOT.connected_z) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(i)) == atomL(lst_positive_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.NOT.connected_z) CYCLE

         if (connected_x.AND.connected_y.AND.connected_z) EXIT
      end do

      do j = 1,natom
         if (atomL(j) == atomL(lst_positive_x%i(i))) then
            write(16,*) atom_name(atom(j)),rxyz(j,1),rxyz(j,2),rxyz(j,3)
         end if
      end do
      CLOSE(UNIT = 16,STATUS = 'KEEP')

CONTAINS

   pure FUNCTION len_3d(r)
      real(wp):: len_3d
      real(wp),intent(in):: r(3)
        len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
   END FUNCTION

END PROGRAM SIMBOX

