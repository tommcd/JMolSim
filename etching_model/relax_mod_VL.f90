MODULE relax_mod
   USE precision_mod, only: wp
   implicit none
   public:: RELAX, RELAX_LIST, RELAX_NMOVE
CONTAINS

   SUBROUTINE RELAX(nr,ifirst,ilast,kbT)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE Keating_energy_mod
      USE repul_force_mod
      USE verlet_list_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,ifirst,ilast
      real(wp),intent(in):: kbT
!     real(wp),parameter:: de = 0.006_wp
      real(wp):: de = 0.08_wp
      integer:: ir,ia,j,k,nfail,ntot
      real(wp):: coprxyz(3),energy1,energy2
      real(wp):: del_energy,dU_kT,erepul
!
      nfail = 0
      ntot = 0
      relax_loop: do ir = 1,NR
      atom_loop: do ia = ifirst,ilast

         if (update(ia)) call new_vlist
         coprxyz(1:3) = rxyz(ia,1:3)

         ! energy of ia before move
         energy1 = energy4(ia,ia)
         do j = 1,ncmax(atom(ia))
            k = proximity(ia,j)
            if (k > 0) then
               energy1 = energy1 + energy4(k,k)
            end if
         end do
         call repul_energy(ia,ia,erepul)
         energy1 = energy1 + erepul

         ! move ia
         do
            rxyz(ia,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
            if (rxyz(ia,3) > 0.0_wp) EXIT
         end do
         rxyz(ia,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         call pbc(rxyz(ia,:))

         ! check if Verlet list should be updated
         if (update(ia)) call NEW_VLIST
         if (update(ia)) then
            write (6, *) 'ERROR: Displacement too large for Verlet '
            STOP
         end if

         ! energy after the move
         energy2 = energy4(ia,ia)
         do j = 1,ncmax(atom(ia))
            k = proximity(ia,j)
            if (k > 0) then
               energy2 = energy2 + energy4(k,k)
            end if
         end do
         call repul_energy(ia,ia,erepul)
         energy2 = energy2 + erepul

         ! Metropolis Monte Carlo decision
         del_energy = energy2 - energy1
         dU_kT = del_energy/kbT
         if (dU_kT <  0.0_wp) GOTO 22
         if (dU_kT > 50.0_wp) GOTO 11
         if (rand() < exp(-dU_kT)) GOTO 22
11       CONTINUE
            ! move is rejected
            rxyz(ia,1:3) = coprxyz(1:3)
            nfail = nfail + 1
22       CONTINUE
         ! move is accepted
         ntot = ntot + 1

      end do atom_loop
      ! write(888, *) energy1,energy2,del_energy
      end do relax_loop
!
      if (ntot == 0) RETURN
      if (real(nfail)/ntot < 0.5_wp) then
         de = de*1.02_wp
      else
         de = de/1.02_wp
      end if
!write(999,'(f0.6,2(i7,1x),f0.6)') de,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX


   SUBROUTINE RELAX_NMOVE(nr,ifirst,ilast,kbT)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE force_keating_mod
      USE repul_force_mod
      USE verlet_list_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,ifirst,ilast
      real(wp),intent(in):: kbT
!     real(wp),parameter:: de = 0.006_wp
      real(wp):: de = 0.01_wp
      integer:: ir,ia,j,k,nfail,ntot
      real(wp):: coprxyz(ifirst:ilast,3)
      real(wp):: ebond,eangle,energy1,energy2
      real(wp):: del_energy,dU_kT,erepul
!
      nfail = 0
      ntot = 0
      relax_loop: do ir = 1,NR

         do ia = ifirst,ilast
            if (update(ia)) then
               call new_vlist
               exit
            end if
         end do
         coprxyz(ifirst:ilast,1:3) = rxyz(ifirst:ilast,1:3)

         ! energy before move
         call energy_keating(ebond,eangle)
         call repul_energy(1,natom,erepul)
         energy1 = ebond + eangle + erepul

         ! move all the atoms
         do ia = ifirst,ilast
            rxyz(ia,1) = coprxyz(ia,1) + (2.0_wp*rand() - 1.0_wp)*de
            rxyz(ia,2) = coprxyz(ia,2) + (2.0_wp*rand() - 1.0_wp)*de
            rxyz(ia,3) = coprxyz(ia,3) + (2.0_wp*rand() - 1.0_wp)*de
         end do
         call pbc(rxyz(ifirst:ilast,:))

         ! check if Verlet list should be updated
         do ia = ifirst,ilast
            if (update(ia)) then
               call NEW_VLIST
               if (update(ia)) then
                  write (6, *) 'ERROR: Displacement too large for Verlet '
                  STOP
               end if
               exit
            end if
         end do

         ! energy after the move
         call energy_keating(ebond,eangle)
         call repul_energy(1,natom,erepul)
         energy2 = ebond + eangle + erepul

         ! Metropolis Monte Carlo decision
         del_energy = energy2 - energy1
         dU_kT = del_energy/kbT
         if (dU_kT <  0.0_wp) GOTO 22
         if (dU_kT > 50.0_wp) GOTO 11
         if (rand() < exp(-dU_kT)) GOTO 22
11       CONTINUE
            ! move is rejected
            rxyz(ifirst:ilast,1:3) = coprxyz(ifirst:ilast,1:3)
            nfail = nfail + 1
22       CONTINUE
         ! move is accepted
         ntot = ntot + 1

      ! write(888, *) energy1,energy2,del_energy
      end do relax_loop
!
      if (ntot == 0) RETURN
      if (real(nfail)/ntot < 0.5_wp) then
         de = de*1.02_wp
      else
         de = de/1.02_wp
      end if
!write(999,'(f0.6,2(i7,1x),f0.6)') de,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX_NMOVE


   SUBROUTINE RELAX_LIST(nr,nl,alist,kbT)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE Keating_energy_mod
      USE repul_force_mod
      USE verlet_list_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,nl,alist(:)
      real(wp),intent(in):: kbT
!     real(wp),parameter:: de = 0.006_wp
      real(wp):: de = 0.08_wp
      integer:: ir,ia,ii,j,k,nfail,ntot
      real(wp):: coprxyz(3),energy1,energy2
      real(wp):: del_energy,dU_kT,erepul
!
      nfail = 0
      ntot = 0
      do ir = 1,NR
      do ii = 1,nl
         ia = alist(ii)
         if (update(ia)) call new_vlist
         coprxyz(1:3) = rxyz(ia,1:3)
         energy1 = energy4(ia,ia)
         do j = 1,ncmax(atom(ia))
            k = proximity(ia,j)
            if (atom(k) == iHydrogen) CYCLE
            if (k > 0) then
               energy1 = energy1 + energy4(k,k)
            end if
         end do
         call repul_energy(ia,ia,erepul)
         energy1 = energy1 + erepul

         do
            rxyz(ia,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
            if (rxyz(ia,3) > 0.0_wp) EXIT
         end do
         rxyz(ia,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         call pbc(rxyz(ia,:))

         if (update(ia)) call NEW_VLIST
         if (update(ia)) then
            write (6, *) 'ERROR: Displacement too large for Verlet '
            STOP
         end if
!
         energy2 = energy4(ia,ia)
         do j = 1,ncmax(atom(ia))
            k = proximity(ia,j)
            if (atom(k) == iHydrogen) CYCLE
            if (k > 0) then
               energy2 = energy2 + energy4(k,k)
            end if
         end do
         call repul_energy(ia,ia,erepul)
         energy2 = energy2 + erepul
!
         del_energy = energy2 - energy1
         dU_kT = del_energy/kbT

         if (dU_kT <  0.0_wp) GOTO 22
         if (dU_kT > 50.0_wp) GOTO 11
         if (rand() < exp(-dU_kT)) GOTO 22
11       CONTINUE
            rxyz(ia,1:3) = coprxyz(1:3)
            nfail = nfail + 1
22       CONTINUE
         ntot = ntot + 1
      end do
!write(888, *) energy1,energy2,del_energy
      end do
      if (real(nfail)/ntot < 0.5_wp) then
         de = de*1.02_wp
      else
         de = de/1.02_wp
      end if
!write(999,'(f0.6,2(i7,1x),f0.6)') de,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX_LIST

END MODULE relax_mod

