
MODULE bond_angle_list_mod
   USE precision_mod
   USE global_vars_mod
   USE bond_angle_types_mod
   USE connectivity_mod, only: proximity
   USE atom_types_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   real(wp),allocatable:: kbond(:),abond(:)
   integer:: nbondtot,nang

CONTAINS

   SUBROUTINE bond_list()
      integer:: ia1,ia2,j
      nbondtot = 0
      if (.not.allocated(nbond)) then
         allocate( ibond(2,natom_max*4),nbond(natom_max) )
         allocate( iang(3,natom_max*6) )
         allocate( kbond(natom_max*4),abond(natom_max*4))
         allocate( kang(natom_max*6),ctheta(natom_max*6))
      end if
      do ia1 = 1,natom
         nbond(ia1) = 0
         do j = 1,ncmax(atom(ia1))
            ia2 = proximity(ia1,j)
            if (ia2 > 0) nbond(ia1) = nbond(ia1) + 1
            if (ia2 > ia1) then
               nbondtot = nbondtot + 1
               ibond(1,nbondtot) = ia1
               ibond(2,nbondtot) = ia2
               call bond_type(ia1,ia2,kbond(nbondtot),abond(nbondtot))
            end if
         end do
      end do
   END SUBROUTINE bond_list


   SUBROUTINE bond_angle_list()
      integer:: ia1,ia2,ia3,ib,ic
      if (.not.allocated(ibond)) then
         allocate( ibond(2,natom_max*4) )
         allocate( iang(3,natom_max*6) )
         allocate( kbond(natom_max*4),abond(natom_max*4))
         allocate( kang(natom_max*6),ctheta(natom_max*6))
      end if
      nbondtot = 0
      nang = 0
      do ia1 = 1,natom
      do ib = 1,ncmax(atom(ia1))
         ia2 = proximity(ia1,ib)
         if (ia2 == 0) cycle
         if (ia1 < ia2) then
            nbondtot = nbondtot + 1
            ibond(1,nbondtot) = ia1
            ibond(2,nbondtot) = ia2
            call bond_type(ia1,ia2,kbond(nbondtot),abond(nbondtot))
            ! print *,ia1,ia2
            ! print *,kbond(nbondtot),abond(nbondtot)
         end if
         do ic = 1,ncmax(atom(ia2))
            ia3 = proximity(ia2,ic)
            if (ia3 == 0) cycle
            if (ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
            if (ia1 < ia3) then
               nang = nang + 1
               iang(1,nang) = ia1
               iang(2,nang) = ia2
               iang(3,nang) = ia3
               call angle_type(ia1,ia2,ia3,kang(nang),ctheta(nang))
            end if
         end do
      end do
      end do
   END SUBROUTINE bond_angle_list


   SUBROUTINE atom_bond_angle(ia,nbnd,ibnd,kbnd,abnd,nan,ian,kan,can)
      integer,intent(in):: ia
      integer,intent(out):: nbnd,ibnd(:,:),nan,ian(:,:)
      real(wp),intent(out):: kbnd(:),abnd(:),kan(:),can(:)
      integer:: L,M,j,k
! the bonds involving atom ia
      nbnd = 0
      do L = 1, ncmax(atom(ia))
         j = proximity(ia,L)
         if (j == 0) cycle
         nbnd = nbnd + 1
         ibnd(1,nbnd) = ia
         ibnd(2,nbnd) = j
         call bond_type(ia,j,kbnd(nbnd),abnd(nbnd))
      end do

! the angles with ia as the centre atom
      nan = 0
      do L = 1, ncmax(atom(ia)) - 1
         j = proximity(ia,L)
         if (j == 0) cycle
         do M = L + 1, ncmax(atom(ia))
            k = proximity(ia,M)
            if (k == 0) cycle
            nan = nan + 1
            ian(1,nan) = j
            ian(2,nan) = ia
            ian(3,nan) = k
            call angle_type(j,ia,k,kan(nan),can(nan))
         end do
      end do

! the angles with ia as an end atom
      do L = 1, ncmax(atom(ia))
         j = proximity(ia,L)
         if (j == 0) cycle
         do M = 1, ncmax(atom(j))
            k = proximity(j,M)
            if (k == ia .or. k == 0) cycle
            nan = nan + 1
            ian(1,nan) = ia
            ian(2,nan) = j
            ian(3,nan) = k
            call angle_type(ia,j,k,kan(nan),can(nan))
         end do
      end do
   END SUBROUTINE atom_bond_angle

END MODULE bond_angle_list_mod

