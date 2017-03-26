MODULE find_list_atoms_mod
   USE precision_mod
   implicit none

CONTAINS

   SUBROUTINE list_atoms_sphere(r0,rad,nat,atomlst)
      USE nlist_mod
      USE coordinates_mod
      real(wp),intent(in):: r0(3),rad
      integer,intent(out):: nat,atomlst(:)
      integer:: ic,nc,j,jj,nl,natmax
      real(wp):: dr(3),r2
! Look for atoms inside the sphere (centre = r0,radius = rad)
      natmax = size(atomlst)
      r2 = rad*rad
      nat = 0
      atomlst = 0
      nl = nlayers_ll(rad)
      if ( nl > nlmax) STOP 'error: list_atoms_sphere, nl > nlmax'
      ic = CELL(r0)
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         j = HOC(nc)
         do while (j /= 0)
            dr = rxyz(j,1:3) - r0
            call pbc(dr(:))
            if ((dr(1)**2 + dr(2)**2 + dr(3)**2) <= r2) then
               nat = nat + 1
               if (nat > natmax) stop 'error: list_atoms_sphere'
               atomlst(nat) = j
            end if
100         j = LL(j)
         end do
      end do cell_loop
   END SUBROUTINE


   SUBROUTINE find_atoms(iat,atomtype,rad,nat,atomlst)
      USE nlist_mod
      USE atom_types_mod
      USE coordinates_mod
      integer,intent(in):: iat,atomtype
      real(wp),intent(in):: rad
      integer,intent(out):: nat,atomlst(:)
      integer:: ic,nc,j,jj,nl
      real(wp):: dr(3),r2
! Look for atoms of type atomtype within a distance rad of ia
      r2 = rad*rad
      nat = 0
      atomlst = 0
      nl = nlayers_ll(rad)
      ic = CELL(rxyz(iat,:))
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
      nc = ncell(jj)
      if (nc == 0) cycle cell_loop
      j = HOC(nc)
      do while (j /= 0)
         if (j == iat) GOTO 100
         if (atom(j) /= atomtype) GOTO 100
         dr(1:3) = rxyz(j,1:3) - rxyz(iat,1:3)
         call pbc(dr(:))
         if ((dr(1)**2 + dr(2)**2 + dr(3)**2) <= r2) then
            nat = nat + 1
            atomlst(nat) = j
         end if
100      j = LL(j)
      end do
      end do cell_loop
   END SUBROUTINE

END MODULE find_list_atoms_mod

