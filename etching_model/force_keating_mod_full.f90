MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none

CONTAINS

   SUBROUTINE energy_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: rl1(3),rl2(3),r1
      real(wp):: rl1l2,dist
!
!
! bond stretching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         dist = ASiO*ASiO
         rl1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
         call pbc(rl1)
         r1 = dot_product(rl1,rl1)
         ebond = ebond + kbond(ib)*(r1 - dist)**2
      end do
      ebond = 0.5_wp*ebond
!
! bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         dist = ASiO*ASiO
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         rl1l2 = dot_product(rl1,rl2)
         eangle = eangle + kang(ia)*(rl1l2 - ctheta(ia)*dist)**2
      end do
      eangle = 0.5_wp*eangle
   END SUBROUTINE energy_keating


   SUBROUTINE force_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: k,rl1(3),r1,ktheta,rl2(3)
      real(wp):: rl1l2,coef,dist
!
! bond stretching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         dist = ASiO*ASiO
         k = KSiO
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         r1 = dot_product(rl1,rl1)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - 2.0_wp*k*(r1 - dist)*rl1(1:3)
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + 2.0_wp*k*(r1 - dist)*rl1(1:3)
         ebond = ebond + 0.5_wp*k*(r1 - dist)**2
      end do
!
! bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         dist = ASiO*ASiO
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         rl1l2 = dot_product(rl1,rl2)
         coef = ktheta*(rl1l2 - ctheta(ia)*dist)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - coef*rl2
         fxyz(ia3,1:3) = fxyz(ia3,1:3) - coef*rl1
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + coef*(rl1 + rl2)
         eangle = eangle + 0.5_wp*ktheta*(rl1l2 - ctheta(ia)*dist)**2
      end do
   END SUBROUTINE force_keating

END MODULE force_keating_mod

