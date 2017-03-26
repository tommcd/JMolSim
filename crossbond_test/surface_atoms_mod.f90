
MODULE SURFACE_ATOMS_MOD
      USE precision_mod, only: wp
      USE coordinates_mod
      USE rand_mod
      implicit none
      logical,allocatable:: surface(:)
      integer,parameter,private:: ntrial = 200
      real(wp):: fsurfa
CONTAINS

   SUBROUTINE SURFACE_ATOMS(ibegin,iend,rprobe)
      USE global_vars, only: natom_max
      USE nlist_mod
      USE atom_types
      USE constants_mod, only: sigma_2
      USE ran_point_sphere,only: ran3sph
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: rprobe
      integer:: ia,it,ic,j,jj,nc
      real(wp):: sx,sy,sz,r0(3)
      real(wp):: dr(3),r(3),rad
      if (.not.allocated(surface)) allocate(surface(natom_max))
      do ia = ibegin,iend
         surface(ia) = .false.
if (atom(ia) == iSilicon) cycle
         rad = rprobe + sigma_2(atom(ia))
         r0(1:3) = rxyz(ia,1:3)
         trial_loop: do it = 1,ntrial
            call ran3sph(sx,sy,sz)
            r(1:3) = r0(1:3) + (/sx,sy,sz/)*rad
            ic = CELL(r)
            CALL NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj = 1,neigh
               nc = ncell(jj)
               if (nc == 0)cycle cell_loop
               j = HOC(nc)
               do while (j /= 0)
if (atom(j) == iSilicon) GOTO 100
                  if (j == ia) GOTO 100
                  dr(1:3) = r(1:3) - rxyz(j,1:3)
                  call pbc(dr)
                  if (dot_product(dr,dr) < (rprobe + sigma_2(atom(j)))**2) then
                     cycle trial_loop
                  end if
100               j = LL(j)
               end do
            end do cell_loop
            surface(ia) = .true.
            exit trial_loop
         end do trial_loop
      end do
      fsurfa = real(count(surface(ibegin:iend)),wp)/(iend-ibegin+1)
   END SUBROUTINE

   SUBROUTINE write_surface_atoms_zbin(iu)
      USE deposition_mod
      USE nlist_mod
      integer,intent(in):: iu
      integer:: iz,ic,nsur,nat,nc,L
      real(wp):: fsur
      do iz = 1,nbin
      nsur = 0
      nat = 0
      CALL Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (surface(L)) then
               nsur = nsur + 1
            end if
            nat = nat + 1
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      fsur = real(nsur,wp)/nat
      write(iu,'(5f12.6)') (iz-0.5_wp)*delz,fsur,fsur/voidage(iz)
      end do
   END SUBROUTINE

END MODULE SURFACE_ATOMS_MOD

