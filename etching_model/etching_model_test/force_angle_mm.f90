module force_keating_mod
   use precision_mod, only: wp
   implicit none
   real(wp):: ebond,eangle
contains

   subroutine force_keating_bond
      use bond_list_mod
      use Keating_parameters
      use coordinates_mod
      integer:: ib,ia1,ia2
      real(wp):: k,dist,rv1(3),r1,df(3)
      ebond = 0.0_wp
!     Force routine for Keating bond streching
      do ib = 1,nbondtot
        k = KSiO
        dist = ASiO
        ia1 = ibond(1,ib)
        ia2 = ibond(2,ib)
        rv1 = rxyz(ia1,1:3)-rxyz(ia2,1:3)
        !call pbc(rv1)
        rv1 = rv1-boxl*anint(rv1*boxli)
        r1 = sqrt(dot_product(rv1,rv1)) 
        df = k*(r1-dist)*rv1(1:3)/r1
        fxyz(ia1,1:3) = fxyz(ia1,1:3) - df
        fxyz(ia2,1:3) = fxyz(ia2,1:3) + df
        ebond = ebond + 0.5_wp*k*(r1-dist)**2
      end do
   end subroutine

   subroutine force_keating_angle
      use bond_list_mod
      use Keating_parameters
      use coordinates_mod
      integer:: ia,i1,i2,i3
      real(wp):: k,cos_theta,cos0,x1,x2,x3,y1,y2,y3,z1,z2,z3
      eangle = 0.0_wp
!     Force routine for Keating bond angle bending
      do ia = 1,nang
         k = kang(ia)
         i1 = iang(1,ia)
         i2 = iang(2,ia)
         i3 = iang(3,ia)

         x2 = rxyz(i2,1)
         y2 = rxyz(i2,2)
         z2 = rxyz(i2,3)

         x1 = rxyz(i1,1)
         y1 = rxyz(i1,2)
         z1 = rxyz(i1,3)

         if((x2 - x1) >  boxl2) x1 = x1 + boxl
         if((x2 - x1) < -boxl2) x1 = x1 - boxl
         if((y2 - y1) >  boxl2) y1 = y1 + boxl
         if((y2 - y1) < -boxl2) y1 = y1 - boxl
         if((z2 - z1) >  boxl2) z1 = z1 + boxl
         if((z2 - z1) < -boxl2) z1 = z1 - boxl

         x3 = rxyz(i3,1)
         y3 = rxyz(i3,2)
         z3 = rxyz(i3,3)

         if((x2 - x3) >  boxl2) x3 = x3 + boxl
         if((x2 - x3) < -boxl2) x3 = x3 - boxl
         if((y2 - y3) >  boxl2) y3 = y3 + boxl
         if((y2 - y3) < -boxl2) y3 = y3 - boxl
         if((z2 - z3) >  boxl2) z3 = z3 + boxl
         if((z2 - z3) < -boxl2) z3 = z3 - boxl

         
         cos_theta = ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                   (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                    Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))
         
         eangle = eangle + 0.5_wp*k*(cos_theta-ctheta(ia))**2
         cos0 = ctheta(ia)
         fxyz(i1,1) = fxyz(i1,1) &
                 -(k*((-x2 + x3)/(Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) -  &
                 ((x1 - x2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))
         fxyz(i1,2) = fxyz(i1,2)  &
                 -(k*((-y2 + y3)/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) - &
                 ((y1 - y2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)*&
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i1,3) = fxyz(i1,3)  &
                 -(k*(-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/  &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 ((-z2 + z3)/(Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) -  &
                 ((z1 - z2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i2,1) = fxyz(i2,1)  &
                 -(k*(((-x2 + x3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5) +  &
                 (-x1 + 2*x2 - x3)/(Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) +  &
                 ((x1 - x2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i2,2) = fxyz(i2,2)  &
                 -(k*(((-y2 + y3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5) +  &
                 (-y1 + 2*y2 - y3)/(Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) +  &
                 ((y1 - y2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i2,3) = fxyz(i2,3)  &
                 -(k*(-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (((-z2 + z3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5) +  &
                 (-z1 + 2*z2 - z3)/(Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)) +  &
                 ((z1 - z2)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**1.5* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i3,1) = fxyz(i3,1)  &
                 -(k*(-(((-x2 + x3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5)) + (x1 - x2)/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i3,2) = fxyz(i3,2)  &
                 -(k*(-(((-y2 + y3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5)) + (y1 - y2)/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))

         fxyz(i3,3) = fxyz(i3,3)  &
                 -(k*(-(((-z2 + z3)*((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3)))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 ((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)**1.5)) + (z1 - z2)/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2)))* &
                 (-cos0 + ((x1 - x2)*(-x2 + x3) + (y1 - y2)*(-y2 + y3) + (z1 - z2)*(-z2 + z3))/ &
                 (Sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)* &
                 Sqrt((-x2 + x3)**2 + (-y2 + y3)**2 + (-z2 + z3)**2))))
      end do
   end subroutine

end module force_keating_mod
