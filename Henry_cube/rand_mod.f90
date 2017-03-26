
MODULE rand_mod
! a random number generator (Numerical Recipies)
! with auxillary subroutines for storing and setting the
! state space
   use precision_mod
   implicit none
   private
   public:: rand,get_rand_state,set_rand_state,read_rand_state,write_rand_state
   integer(i4b),private,save:: ix = -1, iy = -1
   integer(i4b),public,save:: irand_seed = -1
!
   contains

   subroutine set_rand_state(ix0,iy0)
      integer(i4b),intent(in):: ix0,iy0
      real(wp):: tmp
      tmp = rand()
      ix = ix0
      iy = iy0
   end subroutine set_rand_state

   subroutine get_rand_state(ix0,iy0)
      integer(i4b),intent(out):: ix0,iy0
      ix0 = ix
      iy0 = iy
   end subroutine get_rand_state

   subroutine read_rand_state(iu)
      integer(i4b),intent(in):: iu
      real(wp):: tmp
      tmp = rand()
      read(iu,*) ix,iy
      rewind(iu)
   end subroutine read_rand_state

   subroutine write_rand_state(iu)
      integer(i4b),intent(in):: iu
      write(iu,*) ix,iy
      call flush(iu)
   end subroutine write_rand_state

   function rand()
!-----random number generator
      real(wp):: rand
      integer(i4b),parameter::ia=16807,im=2147483647,iq=127773,ir=2836
      integer(i4b),save:: k
      real(dp),save::am
      if(irand_seed <= 0 .or. iy < 0 )then
         am = nearest(1.0_wp,-1.0_wp)/im
         iy = ior(ieor(888889999,abs(irand_seed)),1)
         ix = ieor(777755555,abs(irand_seed))
         irand_seed = abs(irand_seed)+1
      end if
      ix = ieor(ix,ishft(ix,13))
      ix = ieor(ix,ishft(ix,-17))
      ix = ieor(ix,ishft(ix,5))
      k = iy/iq
      iy = ia*(iy-k*iq)-ir*k
      if(iy < 0)iy = iy+im
      rand = am*ior(iand(im,ieor(ix,iy)),1)
   end function rand

! 'minimal standard' RNG
!   function rand0()
!      real(wp):: rand0
!      real(wp),parameter:: im=2147483647.0_wp,ia=16807.0_wp
!      real(wp),parameter:: am=1.0_wp/im
!      irand_seed = mod(ia*mod(ia*irand_seed,im),im)
!      rand0 = am*irand_seed
!   end function

END MODULE rand_mod

