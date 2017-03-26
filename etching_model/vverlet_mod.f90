
MODULE vverlet_mod
    USE precision_mod, only: wp
    implicit none
CONTAINS

! First part of water/salt velocity Verlet algorithm
   SUBROUTINE vverlet_a(dt,massi,r,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: r(:,:),v(:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      r = r + dt*v + dtsq2*f*massi
      v = v + dt2*f*massi
      RETURN
   END SUBROUTINE vverlet_a

! Second part of velocity Verlet
   SUBROUTINE vverlet_b(dt,massi,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: v(:,:)
      v = v + dt*0.5_wp*f*massi
   END SUBROUTINE vverlet_b

END MODULE vverlet_mod

