

program test_force
      integer,parameter:: wp = kind(1.0d0)
      integer:: ia1,ia2,ia3
      real(wp),allocatable:: rxyz(:,:)
      real(wp),allocatable:: fxyz(:,:)
      real(wp):: ktheta
      real(wp):: rl1(3),rl2(3)
      real(wp):: rl1sq,rl2sq,rl1l2,cos_theta,energy1,coef
      real(wp):: coef1,coef2,coef3,energy2,h
!
      allocate(rxyz(3,3),fxyz(3,3))
!
! Force routine for Keating bond angle bending
      h = 0.00000000001_wp
      ktheta = 4.32_wp
      ia1 = 1
      ia2 = 2
      ia3 = 3

      rxyz(1,1) = 0.0_wp
      rxyz(1,2) = 0.161_wp*sqrt(2.0_wp/3.0_wp) + 0.1_wp
      rxyz(1,3) = 0.0_wp

      rxyz(2,1) = 0.161_wp*sqrt(1.0_wp/3.0_wp)
      rxyz(2,2) = 0.0_wp
      rxyz(2,3) = 0.0_wp

      rxyz(3,1) = 0.0_wp
      rxyz(3,2) = -0.161_wp*sqrt(2.0_wp/3.0_wp) - 0.1_wp
      rxyz(3,3) = 0.0_wp

      rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
      rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
!
      rl1sq = sqrt(dot_product(rl1,rl1))
      rl2sq = sqrt(dot_product(rl2,rl2))
      rl1l2 = dot_product(rl1,rl2)
      cos_theta = (rl1l2)/(rl1sq*rl2sq)

      energy1 = 0.5_wp*ktheta*(cos_theta + 1.0_wp/3.0_wp)**2

      coef = ktheta*(cos_theta + 1.0_wp/3.0_wp)
      coef1 = rl1sq**3*rl2sq
      coef2 = rl1sq*rl2sq
      coef3 = rl1sq*rl2sq**3
!
      fxyz(ia2,1:3) = -coef*(-(rl1 + rl2)/coef2 &
                            + rl1l2*rl1/coef1 &
                            + rl1l2*rl2/coef3)
      fxyz(ia1,1:3) = -coef*(rl1/coef2 - rl1l2*rl1/coef1)
      fxyz(ia3,1:3) = -coef*(rl2/coef2 - rl1l2*rl2/coef3)
!-------------------------------------------------------------------


      rxyz(2,1) =  rxyz(2,1) + h




      rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
      rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
!
      rl1sq = sqrt(dot_product(rl1,rl1))
      rl2sq = sqrt(dot_product(rl2,rl2))
      rl1l2 = dot_product(rl1,rl2)
      cos_theta = (rl1l2)/(rl1sq*rl2sq)

      energy2 = 0.5_wp*ktheta*(cos_theta + 1.0_wp/3.0_wp)**2

      write(*,*) (energy2 - energy1)/h, fxyz(ia2,1)
end program test_force


