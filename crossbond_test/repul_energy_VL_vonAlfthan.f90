
MODULE repul_energy_mod
   USE precision_mod, only: wp
   implicit none
   PUBLIC :: repul_energy
CONTAINS

   PURE FUNCTION repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE VERLET_LIST
      real(wp):: repul_energy
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2
      integer:: M,L,jj
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repul_energy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
         if (atom(L) == iHydrogen) CYCLE outer_atom_loop
         atom_loop: do jj = 1, NLIST(L)  ! Verlet-list
            M = LIST(L,jj)
            if (atom(M) == iHydrogen) CYCLE atom_loop
            if (nearest_neighbor2(L,M)) CYCLE atom_loop
            if (M == L) CYCLE atom_loop
            rr(1:3) = rxyz(M,1:3) - rxyz(L,1:3)
            call pbc(rr)
            r2 = rr(1)**2+rr(2)**2+rr(3)**2
            if (r2 <= rden2) repul_energy = repul_energy + 0.5_wp*ff*(r2 - rden2)**2
         end do atom_loop
      end do outer_atom_loop
   END FUNCTION repul_energy

END MODULE repul_energy_mod

