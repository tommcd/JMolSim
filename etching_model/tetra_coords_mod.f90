MODULE tetra_coords_mod
      USE precision_mod
      implicit none
CONTAINS

   SUBROUTINE tetra_coords(bondl,tetra)
      real(wp),intent(in):: bondl
      real(wp),intent(out):: tetra(3,4)
!
! returns a set of coordinates for the 4 outer atoms of a
! tetrahedral molecule (with the central atom at the origin).
! Atom 1 is on the z axis below the origin and atoms 2,3, & 4
! are on a plane above the origin.
!
      real(wp):: zz,sth,cth,rxy,xx,yy
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,2) = (/ xx, yy, zz /)
      tetra(1:3,3) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,4) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

   SUBROUTINE tetra_coords5(bondl,tetra)
      real(wp),intent(in):: bondl
      real(wp),intent(out):: tetra(3,5)
!
! returns a set of coordinates for a tetrahedral molecule
! with the central atom at the origin (atom 1), one atom
! on the z axis below it (atom 2), and atoms 3,4, and 5 on
! a plane above the origin.
!
      real(wp):: zz,sth,cth,rxy,xx,yy
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

   SUBROUTINE SiOH4_coords(bondl,bondl_OH,angle_SiOH,tetra)
      USE rotate_axis_mod
      USE rand_mod
      real(wp),intent(in):: bondl,bondl_OH,angle_SiOH
      real(wp),intent(out):: tetra(3,9)
!
! Returns a set of coordinates for a tetrahedral molecule
! with the central atom at the origin (atom 1), one atom
! on the z axis below it (atom 2), and atoms 3,4, and 5 on
! a plane above the origin.
! Atoms 6,7,8 and 9 are attached to atoms 2,3,4 and 5
! respectively
!
      real(wp):: zz,sth,cth,rxy,xx,yy,rr(3)
      real(wp):: aa(3,3),cang,sang,pi,phi
      integer:: j
      pi = 4.0_wp*(atan(1.0_wp))
      cang = cos(angle_SiOH)
      sang = sin(angle_SiOH)
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
! attach atoms 6-9 to 2-5
      do j = 2,5
         rr = tetra(1:3,j)/bondl
         call zaxis2vect(rr,aa)
         phi = rand()*2.0_wp*pi
         zz = -bondl_oh*cang
         xx = bondl_oh*sang*cos(phi)
         yy = bondl_oh*sang*sin(phi)
         tetra(1:3,j + 4) = tetra(1:3,j) + matmul(aa,(/ xx,yy,zz /))
      end do
   END SUBROUTINE

END MODULE

