
MODULE repul_energy_mod
   USE precision_mod, only: wp
   implicit none
   PUBLIC :: repul_energy
CONTAINS

   PURE FUNCTION repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      real(wp):: repul_energy
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2
      integer:: M,L,ic,nc,inn
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repul_energy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      if (atom(L) == iHydrogen) CYCLE outer_atom_loop
      ic = CELL(rxyz(L,:)) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if (atom(M) == iHydrogen) GOTO 100
            if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            rr(1:3) = rxyz(M,1:3) - rxyz(L,1:3)
            call pbc(rr)
            r2 = rr(1)**2+rr(2)**2+rr(3)**2
            if (r2 <= rden2) then
               repul_energy = repul_energy + 0.5_wp*ff*(r2 - rden2)**2
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
   END FUNCTION repul_energy

END MODULE repul_energy_mod

