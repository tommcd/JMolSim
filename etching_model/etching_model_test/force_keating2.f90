
module force_keating_mod
   use precision_mod, only: wp
   implicit none
   real(wp):: energy
contains
   subroutine force_keating
      use bond_list_mod
      use Keating_parameters
      use coordinates_mod
      integer:: ib,ia,ia1,ia2,ia3,i_nonb,j
      real(wp):: k,dist,rl1(3),rl1sq,ktheta,rl2(3),rl2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
      fxyz = 0.0_wp
      energy = 0.0_wp

!     Force routine for Keating bond streching

      do ib = 1,nbondtot
        k = KSiO
        dist = ASiO
        ia1 = ibond(1,ib)
        ia2 = ibond(2,ib)
        rl1 = rxyz(ia1,1:3)-rxyz(ia2,1:3)
        !call pbc(rl1)
        rl1 = rl1-boxl*anint(rl1*boxli)

        rl1sq = sqrt(dot_product(rl1,rl1)) 

        fxyz(ia1,1:3) = fxyz(ia1,1:3)-k*(rl1sq-dist)*rl1(1:3)/rl1sq
        fxyz(ia2,1:3) = fxyz(ia2,1:3)+k*(rl1sq-dist)*rl1(1:3)/rl1sq

        energy = energy + 0.5_wp*k*(rl1sq-dist)**2
        
      end do

!     Force routine for Keating bond angle bending

      do ia = 1,nang
        ktheta = kang(ia)
        ia1 = iang(1,ia)
        ia2 = iang(2,ia)
        ia3 = iang(3,ia)
        rl1 = rxyz(ia1,1:3)-rxyz(ia2,1:3)
        rl2 = rxyz(ia3,1:3)-rxyz(ia2,1:3)
        !call pbc(rl1)
        !call pbc(rl2)
        rl1 = rl1-boxl*anint(rl1*boxli)
        rl2 = rl2-boxl*anint(rl2*boxli)

        rl1sq = sqrt(dot_product(rl1,rl1)) 
        rl2sq = sqrt(dot_product(rl2,rl2))
        rl1l2 = dot_product(rl1,rl2)
        cos_theta = (rl1l2)/(rl1sq*rl2sq)
        coef = ktheta*(cos_theta-ctheta(ia))
        coef1 = rl1sq**3*rl2sq
        coef2 = rl1sq*rl2sq 
        coef3 = rl1sq*rl2sq**3 

        fxyz(ia2,1:3) = fxyz(ia2,1:3)-coef*(-(rl1+rl2)/coef2 &
                                        + rl1l2*rl1/coef1 &
                                        + rl1l2*rl2/coef3)
        fxyz(ia1,1:3) = fxyz(ia1,1:3)-coef*(rl2/coef2-rl1l2*rl1/coef1)
        fxyz(ia3,1:3) = fxyz(ia3,1:3)-coef*(rl1/coef2-rl1l2*rl2/coef3)
        energy = energy + 0.5_wp*ktheta*(cos_theta-ctheta(ia))**2
        
      end do

   end subroutine
end module force_keating_mod
