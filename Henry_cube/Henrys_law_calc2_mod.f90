MODULE Henrys_law_calc_mod
    USE precision_mod, only: wp
    USE global_vars
    implicit none
    PUBLIC:: Henrys_law_calc
CONTAINS

   SUBROUTINE Henrys_law_calc(ntrial,xl,xu,yl,yu,zl,zu,pr,fover,Uljel,Khenry,facc)
      USE coordinates_mod
      USE atom_types
!     USE constants_mod, only: sigma_2
      USE seaton_mod, only: sigi2
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE global_vars, only: kboltzT
      USE lj_el_mod
      USE probe_mol_mod
      integer,intent(in):: ntrial
      real(wp),intent(in):: xl,xu,yl,yu,zl,zu
      type(probe_mol),intent(in):: pr
      real(wp),intent(in):: fover
      real(wp),intent(out):: Uljel,Khenry,facc
!     real,parameter:: expon_max = log(huge(1.0_wp))*0.5_wp  ! not standard f95
      real:: expon_max
      type(probe_mol):: pr1
      real(wp):: aa(3,3),q(4),x,y,z,rnaccept,dUlj,dUel,dKHenry,expon
      integer:: it,nat,i
      logical:: rotate
      expon_max = log(huge(1.0_wp))*0.5_wp
      nat = pr%n
      pr1 = pr
      rotate = ((size(pr%r,2) > 1) .and. nat > 1)
      rnaccept = 0
      Uljel = 0.0_wp
      KHenry = 0.0_wp
      insert_loop: do it = 1,ntrial
         x = xl + rand()*(xu - xl)
         y = yl + rand()*(yu - yl)
         z = zl + rand()*(zu - zl)
         if (rotate) then
            !call random_rotation(aa)
            call ran4sph(q)
            call quat2mat(q,aa)
            pr1%r = matmul(aa,pr%r)
         else
            pr1%r = pr%r
         end if
         pr1%r(1,:) = pr1%r(1,:) + x
         pr1%r(2,:) = pr1%r(2,:) + y
         pr1%r(3,:) = pr1%r(3,:) + z
!!       pr1%r(1:2,:) = pr1%r(1:2,:) - anint(pr1%r(1:2,:)*boxli)*boxl
         do i = 1,nat
            call pbc(pr1%r(1:3,i))
         end do
! first check if there is a large overlap
         if(overlap(fover)) cycle insert_loop
         rnaccept = rnaccept + 1.0_wp
         call ENERGY_LJ_EL2(pr1,2,dUlj,dUel)
!print *,dUlj,dUel
         Uljel = Uljel + dUlj + dUel
         expon = -(dUlj+dUel)/kboltzT
         expon = min(expon,expon_max)
         dKHenry = exp(expon)
         KHenry = KHenry + dKHenry
      end do insert_loop
      Uljel = Uljel/real(ntrial,wp)
      KHenry = KHenry/real(ntrial,wp)
      facc = rnaccept/real(ntrial,wp)
      RETURN
!
   contains
!
      pure function overlap(fr)
         real(wp),intent(in):: fr
         logical:: overlap
         integer:: k,jj,j,ic,nc,ncell(125),neigh
         real(wp):: dr(3)
         overlap = .false.
         probe_atom_loop: do k = 1,nat
         ic = CELL( pr1%r(1:3,k) )
         call NEIGCELL(ic,1,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0) cycle cell_loop
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               if (atom(j) /= iSilicon) then
               dr(1:3) = pr1%r(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               if (dot_product(dr,dr) < ( fr*( pr1%rad(k) + sigi2(atom(j)) ) )**2) then
                  overlap = .true.
                  return
               end if
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
         end do probe_atom_loop
      end function overlap
   END SUBROUTINE Henrys_law_calc

END MODULE Henrys_law_calc_mod

