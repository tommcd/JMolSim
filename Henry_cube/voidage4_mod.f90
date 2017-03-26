
MODULE voidage_mod
   USE precision_mod
   implicit none
CONTAINS

   FUNCTION voidage_calc(ntrial,zlower,zupper,pr)
      USE coordinates_mod
      USE atom_types, only: atom,iSilicon
      USE constants_mod, only: sigma_2
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE probe_mol_mod
      real(wp):: voidage_calc
      integer,intent(in):: ntrial
      real(wp),intent(in):: zlower,zupper
      type(probe_mol),intent(in):: pr
      real(wp),allocatable:: probe(:,:)
      real(wp):: aa(3,3),q(4)
      real(wp):: x,y,z,dr(3),rnaccept
      integer:: i,k,j,ic,jj,nc,it,nat
      logical:: rotate
      nat = pr%n
      allocate(probe(3,nat))
      rotate = ((size(probe,2) > 1) .and. nat > 1)
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         z = zlower + rand()*(zupper - zlower)
         if (rotate) then
            !call random_rotation(aa)
            call ran4sph(q)
            call quat2mat(q,aa)
            probe = matmul(aa,pr%r)
         else
            probe = pr%r
         end if
         probe(1,1:nat) = probe(1,1:nat) + x
         probe(2,1:nat) = probe(2,1:nat) + y
         probe(3,1:nat) = probe(3,1:nat) + z
         do i=1,nat
            call pbc(probe(1:3,i))
         end do
! check for overlap
         probe_atom_loop: do k = 1,nat
            ic = CELL( probe(1:3,k) )
            call NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj = 1,neigh
               nc = ncell(jj)
               if (nc == 0)cycle cell_loop
               j = HOC(nc)
               cell_atom_loop: do while (j /= 0)
                  if (atom(j) /= iSilicon) then
                  dr(1:3) = probe(1:3,k) - rxyz(j,1:3)
                  call pbc(dr)
                  if (dot_product(dr,dr) < (pr%rad(k) + sigma_2(atom(j)))**2) then
                     cycle insert_loop
                  end if
                  end if
                  j = LL(j)
               end do cell_atom_loop
            end do cell_loop
         end do probe_atom_loop
         rnaccept = rnaccept + 1.0_wp
      end do insert_loop
      voidage_calc = rnaccept/real(ntrial,wp)
   END FUNCTION voidage_calc

   FUNCTION voidage_calc1(ntrial,zlower,zupper,rprobe)
      USE coordinates_mod
      USE atom_types, only: atom,iSilicon
      USE constants_mod, only: sigma_2
      USE rand_mod
      USE nlist_mod
      real(wp):: voidage_calc1
      integer,intent(in):: ntrial
      real(wp),intent(in):: zlower,zupper,rprobe
      real(wp):: r(3),dr(3),rnaccept
      integer:: j,ic,jj,nc,it
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         r(1) = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         r(2) = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         r(3) = zlower + rand()*(zupper - zlower)
         ic = CELL(r)
         call NEIGCELL(ic,1,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0)cycle cell_loop
            j = HOC(nc)
            do while (j /= 0)
               if (atom(j) /= iSilicon) then
                  dr(1:3) = r(1:3) - rxyz(j,1:3)
                  call pbc(dr)
                  if (dot_product(dr,dr) < (rprobe + sigma_2(atom(j)))**2) then
                     cycle insert_loop
                  end if
               end if
               j = LL(j)
            end do
         end do cell_loop
         rnaccept = rnaccept + 1.0_wp
      end do insert_loop
      voidage_calc1 = rnaccept/real(ntrial,wp)
   END FUNCTION voidage_calc1

END MODULE voidage_mod

