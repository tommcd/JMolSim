
MODULE bond_list_mod
   use precision_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   integer:: nbondtot,nang
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types
      USE Keating_parameters
      integer:: ia1,ia2,j
      nbondtot = 0
      if (.not.allocated(nbond)) then
         allocate( ibond(2,natom_max*4),nbond(natom_max) )
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
            end if
         end do
      end do
   END SUBROUTINE bond_list
   
   subroutine set_bond_angle_lists()
      USE global_vars, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types
      USE Keating_parameters
      integer:: i,j,ia1,ia2,ia3,ib,ic
!
      if (.not.allocated(ibond)) then
         allocate( ibond(2,natom_max*4) )
         allocate( iang(3,natom_max*6) )
      end if
      allocate(kang(6*natom_max),ctheta(6*natom_max))
      nbondtot = 0
      nang = 0
      do ia1 = 1,natom
      do ib = 1,ncmax(atom(ia1))
         ia2 = proximity(ia1,ib)
         if(ia2 == 0)cycle
         if(ia1 < ia2)then
            nbondtot = nbondtot+1
            ibond(1,nbondtot) = ia1
            ibond(2,nbondtot) = ia2
         end if
         do ic = 1,ncmax(atom(ia2))
            ia3 = proximity(ia2,ic)
            if(ia3 == 0)cycle
            if(ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
            if(ia1 < ia3)then
               nang = nang+1
               iang(1,nang) = ia1
               iang(2,nang) = ia2               
               iang(3,nang) = ia3
            end if
         end do
      end do
      end do
      do i = 1, nang
         if (atom(iang(2,i)) == iSilicon) then
            kang(i)=KOSiO
            ctheta(i)=-1.0_wp/3.0_wp
         else
            kang(i)=KSiOSi
            ctheta(i)=-1.0_wp
         end if
      end do     
!

   end subroutine set_bond_angle_lists
END MODULE bond_list_mod

