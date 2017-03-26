
module vverlet_mod
    use precision_mod, only: wp
    implicit none
contains

! First part of water/salt velocity Verlet algorithm
   subroutine vverlet_a(dt,massi,r,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: r(:,:),v(:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      r = r + dt*v + dtsq2*f*massi
      v = v + dt2*f*massi
      return
   end subroutine vverlet_a

! Second part of velocity Verlet
   subroutine vverlet_b(dt,massi,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: v(:,:)
      v = v + dt*0.5_wp*f*massi
   end subroutine vverlet_b

end module vverlet_mod

