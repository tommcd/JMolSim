
MODULE lj_el_mod
    USE precision_mod, only: wp
    USE global_vars
    implicit none
    private
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: uljcut(:,:),uljsf(:,:)
    real(wp),parameter:: rcut = 12.0_wp/angstrom
    real(wp),parameter:: rcut2 = rcut**2
    PUBLIC:: LJ_INIT,ENERGY_LJ_EL,ENERGY_LJ_EL2
CONTAINS

   SUBROUTINE LJ_INIT
      USE seaton_mod, only: epsi,sigi,ntyplj
      real(wp):: eij,rstij,rcuti,rcuti6
      integer:: i,j
      rcuti = (1.0_wp/rcut)
      rcuti6 = rcuti**6
      if(.not.allocated(aij)) allocate( aij(0:ntyplj,0:ntyplj) )
      if(.not.allocated(bij)) allocate( bij(0:ntyplj,0:ntyplj) )
      if(.not.allocated(uljcut)) allocate( uljcut(0:ntyplj,0:ntyplj) )
      if(.not.allocated(uljsf))  allocate(  uljsf(0:ntyplj,0:ntyplj) )
      do j = 0,ntyplj    ! lj aij,bij parameters & shifted force terms
      do i = 0,ntyplj
        eij = sqrt(epsi(i)*epsi(j))
        rstij = (sigi(i)+sigi(j))*0.5
        aij(i,j) = 4.0*(rstij**12)*eij
        bij(i,j) = 4.0*(rstij**6)*eij
        uljcut(i,j) = rcuti6*(aij(i,j)*rcuti6-bij(i,j))
        uljsf(i,j) = -rcuti6*(12.0*aij(i,j)*rcuti6-6.0*bij(i,j))*rcuti
!        write(*,'(2i6,4g18.9)')i,j,aij(i,j),bij(i,j),uljcut(i,j),uljsf(i,j)
      end do
      end do
   END SUBROUTINE

   PURE SUBROUTINE ENERGY_LJ_EL(pr,sysfirst,syslast,Ulj,Uel)
      USE coordinates_mod
      USE atom_types, only: atom
      USE charges_mod, only: charge
      USE probe_mol_mod
      type(probe_mol),intent(in):: pr
      integer,intent(in):: sysfirst,syslast
      real(wp),intent(out):: Ulj,Uel
      real(wp):: dr(3),r2,r1,rr6,dele,qL
      integer:: L,M,iaL
      Ulj = 0.0_wp
      Uel = 0.0_wp
      outer_atom_loop: do L = 1,pr%n
      iaL = pr%atom(L)
      qL = pr%q(L)
      atom_loop: do M = sysfirst,syslast
         dr(1:3) = rxyz(M,1:3) - pr%r(1:3,L)
         call pbc(dr)
         r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
         r1 = sqrt(r2)
         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(iaL,atom(M))*rr6-bij(iaL,atom(M))) &
                 - uljcut(iaL,atom(M)) - uljsf(iaL,atom(M))*(r1-rcut)
            Ulj = Ulj + dele
         end if
         Uel = Uel + qL*charge(M)/r1
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE

   PURE SUBROUTINE ENERGY_LJ_EL2(pr,nl,Ulj,Uel)
      USE coordinates_mod
      USE atom_types, only: atom
      USE charges_mod, only: charge
      USE probe_mol_mod
      USE nlist_mod,only: cell,neigcell,hoc,ll
      type(probe_mol),intent(in):: pr
      integer,intent(in):: nl
      real(wp),intent(out):: Ulj,Uel
      real(wp):: dr(3),r2,r1,rr6,dele,qL
      integer:: L,M,iaL
      integer:: jj,ic,nc,ncell(125),neigh
      Ulj = 0.0_wp
      Uel = 0.0_wp
      outer_atom_loop: do L = 1,pr%n
      iaL = pr%atom(L)
      qL = pr%q(L)
      ic = CELL( pr%r(1:3,L) )
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         cell_atom_loop: do while (M /= 0)
            dr(1:3) = rxyz(M,1:3) - pr%r(1:3,L)
            call pbc(dr)
            r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
            r1 = sqrt(r2)
            if (r2 < rcut2) then
               rr6 = 1.0_wp/(r2**3)
               dele = rr6*(aij(iaL,atom(M))*rr6-bij(iaL,atom(M))) &
                    - uljcut(iaL,atom(M)) - uljsf(iaL,atom(M))*(r1-rcut)
               Ulj = Ulj + dele
            end if
            Uel = Uel + qL*charge(M)/r1
            M = LL(M)
         end do cell_atom_loop
      end do cell_loop
      end do outer_atom_loop
   end subroutine


   PURE SUBROUTINE EnergyLJEL(atomfirst,atomlast,sysfirst,syslast,Uljel)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE atom_types, only: atom
      USE charges_mod, only: charge
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: Uljel
      real(wp):: dr(3),r2,r1,rr6,dele
      integer:: L,M
      Uljel = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         if (nearest_neighbor2(M,L)) cycle atom_loop
         dr(1:3) = rxyz(M,1:3) - rxyz(L,1:3)
         call pbc(dr)
         r2 = dr(1)**2+dr(2)**2+dr(3)**2
         r1 = sqrt(r2)
         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(atom(L),atom(M))*rr6-bij(atom(L),atom(M))) &
                 - uljcut(atom(L),atom(M)) - uljsf(atom(L),atom(M))*(r1-rcut)
            Uljel = Uljel + dele
         end if
         Uljel = Uljel + charge(L)*charge(M)/r1
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE

END MODULE lj_el_mod

