
module rk_mod
    use precision_mod, only: wp
    implicit none
    real(wp),save,allocatable:: k1_rk(:,:),k2_rk(:,:)
contains

! First part of water/salt velocity Verlet algorithm
   subroutine rk_a(dt,massi,r,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: r(:,:),v(:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      k1_rk = dtsq2*f*massi
      r = r + (2.0_wp/3.0_wp)*dt*v + (4.0_wp/9.0_wp)*k1_rk
      return
   end subroutine

! Second part of velocity Verlet
   subroutine rk_b(dt,massi,r,rold,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:),rold(:,:)
      real(wp),intent(inout):: v(:,:)
      real(wp),intent(out):: r(:,:)
      k2_rk = dt*dt*0.5_wp*f*massi
      r = rold + v*dt + 0.5_wp*(k1_rk+k2_rk)
      v = v + 0.5_wp*(k1_rk+3.0_wp*k2_rk)/dt
   end subroutine

end module rk_mod

