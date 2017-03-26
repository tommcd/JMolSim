
MODULE relax_mod
   USE precision_mod, only: wp
   implicit none
   public:: RELAX
CONTAINS

   SUBROUTINE RELAX(nr,ifirst,ilast)
      USE coordinates_mod
      USE atom_types
      USE connectivity_mod
      USE Keating
      USE repul_energy_mod
      USE verlet_list
      USE rand_mod
      USE global_vars, only: etot
      implicit none
      integer,intent(in):: nr,ifirst,ilast
      real(wp),parameter:: de = 0.006_wp
!      real(wp):: de = 0.08_wp
      integer:: ir,i,j,k,nfail
      real(wp):: coprxyz(3),energy1,energy2,del_energy
      real(wp):: EEE
!

print*, 'etot = ', etot
      nfail = 0

      do ir = 1,NR
print*, 'ir = ', ir

      do i = ifirst,ilast
         if (update(i)) then
            call new_vlist
            print*, 'update = ', update(i)
         end if 
         coprxyz(1:3) = rxyz(i,1:3)
         energy1 = energy4(i,i)

         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (atom(k) == iHydrogen) CYCLE
            if (k > 0) then
               energy1 = energy1 + energy4(k,k) 
            end if
         end do

         energy1 = energy1 + repul_energy(i,i)

!print*, 'rxyz(',i,',:)[before]     = ', rxyz(i,:)
         rxyz(i,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(i,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(i,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
!print*, 'rxyz(',i,',:)[after]      = ', rxyz(i,:)
         call pbc(rxyz(i,:))
         
         if (update(i)) call NEW_VLIST
         if (update(i)) then
            write (6, *) 'ERROR: Displacement too large for Verlet '
            STOP
         end if
!
         energy2 = energy4(i,i)
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (atom(k) == iHydrogen) CYCLE
            if (k > 0) then
               energy2 = energy2 + energy4(k,k)
            end if
         end do
         energy2 = energy2 + repul_energy(i,i)
!
         del_energy = energy2 - energy1
         EEE = del_energy/Etot
         if (EEE <  0.0_wp) GOTO 22
         if (EEE > 50.0_wp) GOTO 11
         if (rand() < exp(-EEE)) GOTO 22
11       CONTINUE
            rxyz(i,1:3) = coprxyz(1:3)
            nfail = nfail + 1
22       CONTINUE
      end do

      end do

!print*, 'nfail = ', nfail
!print*, 'NR*(ilast-ifirst+1) = ', NR*(ilast-ifirst+1)
!print*, 'acceptance ratio = ', real(nfail)/(NR*(ilast-ifirst+1))
!print*, ''

!     print '(f0.6,2(i7,1x),f0.6)',de,NR*(ilast-ifirst+1),nfail,real(nfail)/(NR*(ilast-ifirst+1))
!     if (real(nfail)/(NR*(ilast-ifirst+1)) < 0.5) then
!        de = de*1.01
!     else
!        de = de/1.01
!     end if
   END SUBROUTINE RELAX

END MODULE relax_mod

