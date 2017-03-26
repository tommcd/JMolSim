
MODULE set_bond_angle_lists_mod
   implicit none
!
CONTAINS

   SUBROUTINE set_bond_angle_lists()
      integer:: i,j,ia1,ia2,ia3,ib,ic
!
      nbond = 0
      nang = 0
      do ia1 = 1,natom
      do ib = 1,ncmax(atom(ia1))
         ia2 = proximity(ib,ia1)
         if (ia2 == 0) cycle
         if (ia1 < ia2) then
            nbond = nbond + 1
            ibond(1,nbond) = ia1
            ibond(2,nbond) = ia2
         end if
         do ic = 1,ncmax(atom(ia2))
            ia3 = proximity(ic,ia2)
            if (ia3 == 0) cycle
            if (ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
            if (ia1 < ia3) then
               nang = nang + 1
               iang(1,nang) = ia1
               iang(2,nang) = ia2
               iang(3,nang) = ia3
            end if
         end do
      end do
      end do
!
   END SUBROUTINE set_bond_angle_lists
END MODULE set_bond_angle_lists_mod

