
MODULE ran_point_sphere
   USE precision_mod
   USE rand_mod, only: rand
   implicit none
CONTAINS

   SUBROUTINE ran_point_in_sphere(qx,qy,qz)
      real(wp),intent(out):: qx,qy,qz
      real(wp)::s1
      do
         qx = (2.0_wp*rand()-1.0_wp)
         qy = (2.0_wp*rand()-1.0_wp)
         qz = (2.0_wp*rand()-1.0_wp)
         s1=qx*qx+qy*qy+qz*qz
         if (s1 <= 1.0_wp)exit
      end do
      return
   END SUBROUTINE

   SUBROUTINE ran3sph(x,y,z)
!
! random vector on the surface of a sphere
! Marsaglia,G., Ann.math.Stat.,43,645-6,(1972)
!
      real(wp),intent(out):: x,y,z
      real(wp):: s,tmp

1     continue
         x = 2.0_wp*rand() - 1.0_wp
         y = 2.0_wp*rand() - 1.0_wp
         s = x*x + y*y
      if (s >= 1.0_wp) GOTO 1
      z = 2.0_wp*s - 1.0_wp
      tmp = 2.0_wp*sqrt(1.0_wp - s)
      x = x*tmp
      y = y*tmp
      return
   END SUBROUTINE

   SUBROUTINE ran3sph2(x,y,z)
      real(wp),parameter::twopi=6.283185307179586476925286766559005768394_wp
      real(wp),intent(out)::x,y,z
      real(wp)::s,tmp
      z = 2.0_wp*rand()-1.0_wp
      tmp = rand()*twopi
      s = sqrt(1.0_wp-z*z)
      x = s*cos(tmp)
      y = s*sin(tmp)
      return
   END SUBROUTINE

   SUBROUTINE ran_point_on_sphere(qx,qy,qz)
      real(wp),intent(out):: qx,qy,qz
      real(wp)::s1,r1
      do
        qx = (2.0_wp*rand()-1.0_wp)
        qy = (2.0_wp*rand()-1.0_wp)
        qz = (2.0_wp*rand()-1.0_wp)
        s1=qx*qx+qy*qy+qz*qz
        if (s1 <= 1.0_wp)exit
      end do
      r1 = 1.0_wp/sqrt(s1)
      qx = qx*r1
      qx = qx*r1
      qz = qx*r1
      return
   END SUBROUTINE

   SUBROUTINE ran4sph(q)
! random point on surface of 4-Sphere (a Quaternion)
      real(wp),intent(out)::q(4)
      real(wp)::s1,s2
      do
        q(1) = (2.0_wp*rand()-1.0_wp)
        q(2) = (2.0_wp*rand()-1.0_wp)
        s1 = q(1)*q(1)+q(2)*q(2)
        if (s1 <= 1.0_wp)exit
      end do
      do
        q(3) = (2.0_wp*rand()-1.0_wp)
        q(4) = (2.0_wp*rand()-1.0_wp)
        s2 = q(3)*q(3) + q(4)*q(4)
        if (s2 <= 1.0_wp)exit
      end do
      q(3) = q(3)*sqrt((1.0_wp-s1)/s2)
      q(4) = q(4)*sqrt((1.0_wp-s1)/s2)
      return
   END SUBROUTINE

END MODULE

