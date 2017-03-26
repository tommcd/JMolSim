MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none
CONTAINS

   SUBROUTINE energy_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: rl1(3),r1sq,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         rl1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         ebond = ebond + kbond(ib)*(r1sq - abond(ib))**2
      end do
      ebond = 0.5_wp*ebond
!
!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         eangle = eangle + kang(ia)*(cos_theta - ctheta(ia))**2
      end do
      eangle = 0.5_wp*eangle
   END SUBROUTINE energy_keating


   SUBROUTINE force_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: k,dist,rl1(3),r1sq,ktheta,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         k = KSiO
         dist = ASiO
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - k*(r1sq - dist)*rl1(1:3)/r1sq
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + k*(r1sq - dist)*rl1(1:3)/r1sq
         ebond = ebond + 0.5_wp*k*(r1sq - dist)**2
      end do

!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         coef = ktheta*(cos_theta - ctheta(ia))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
         fxyz(ia2,1:3) = fxyz(ia2,1:3) - coef*(-(rl1 + rl2)/coef2 &
                                            + rl1l2*rl1/coef1 &
                                            + rl1l2*rl2/coef3)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - coef*(rl2/coef2 - rl1l2*rl1/coef1)
         fxyz(ia3,1:3) = fxyz(ia3,1:3) - coef*(rl1/coef2 - rl1l2*rl2/coef3)
         eangle = eangle + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE force_keating

END MODULE force_keating_mod

