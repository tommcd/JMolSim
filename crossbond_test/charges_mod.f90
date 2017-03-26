
MODULE charges_mod
   USE precision_mod
   USE global_vars, only: nseed
   USE atom_types
   USE seaton_mod, only: qi => q_seaton
   implicit none
   PRIVATE
   real(wp),allocatable:: charge(:)
   PUBLIC:: charge, assign_charge, qsi
!
CONTAINS

   pure function qsi(n)
      real(wp):: qsi
      integer,intent(in):: n
      qsi = (n-4)*qi(iOxygen)*0.5_wp - n*(qi(iOxygenH) + qi(iHydrogen))
   end function

   SUBROUTINE assign_charge(ifirst,ilast)
      USE connectivity_mod, only: noh
      USE atom_types
      integer,intent(in):: ifirst,ilast
      real(wp):: qo,qh,qoh
      integer:: i
      qoh= qi(iOxygenH)
      qh = qi(iHydrogen)
      qo = qi(iOxygen)
      do i = ifirst,nseed
         select case(atom(i))
         case(iOxygen)
            charge(i) = qo
            charge(i-nseed) = -qo*0.5_wp
         case(iOxygenH)
            charge(i) = qoh
            charge(i-nseed) = -(qoh + qh)
         case default
            stop 'assign_charge: wrong seed type'
         end select
      end do
      do i = nseed+1,ilast
         select case(atom(i))
         case(iSilicon)
            charge(i) = qsi(noh(i))
         case(iOxygen)
            charge(i) = qo
         case(iOxygenH)
            charge(i) = qoh
         case(iHydrogen)
            charge(i) = qh
         case default
            stop 'assign_charge: unknown atom type'
         end select
      end do
   END SUBROUTINE

END MODULE charges_mod

