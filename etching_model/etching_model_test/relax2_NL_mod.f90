
MODULE relax_mod
   use precision_mod, only: wp
   implicit none
   real(wp):: Etot
   public:: RELAX
CONTAINS

   SUBROUTINE RELAX(nr,ifirst,ilast,delt,alist,na,quench,scale_repul)
      USE BLOCK1
      USE Keating
      USE repul_energy_mod
      use nlist_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,ifirst,ilast,alist(:),na
      logical,intent(in):: quench,scale_repul
      real(wp),intent(in):: delt
      real(wp),PARAMETER:: de = 0.006_wp
!     real(wp):: de = 0.08_wp
      integer:: ir,i,j,k,nfail
      real(wp):: coprxyz(3),energy1,energy2,del_energy
      real(wp):: EEE,fscale0,fscale,E0
      logical:: special_atom,any_special_atoms,update
!
      nfail = 0
      any_special_atoms = ( alist(1)/=0 )
      DO ir = 1,NR

      if(scale_repul) fscale0 = real(ir-1,wp)/real(NR-1,wp)
      if(quench) E0 = Etot + (NR-ir)*delt/(NR-1)

      DO i = ifirst,ilast

         if (any_special_atoms)then
            special_atom = any(alist(1:na) == i)
         else
            special_atom = .false.
         end if
         if(special_atom .and. scale_repul)then
            fscale = fscale0
         else
            fscale = 1.0_wp
         end if

         coprxyz(1:3) = rxyz(i,1:3)
         energy1 = energy4(i,i)
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if(k > 0)then
               energy1 = energy1 + energy4(k,k)
            end if
         end do
         energy1 = energy1 + fscale*repul_energy(i,i)

         DO
            rxyz(i,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
            IF(rxyz(i,3) > 0.0_wp) EXIT
         END DO
         rxyz(i,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(i,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(i,1:2) = rxyz(i,1:2) - boxl*ANINT(rxyz(i,1:2)*boxli)

         update = .false.
         IF (CELL(rxyz(i,1:3)) /= CELL(coprxyz))then
            call NEW_NLIST
            update = .true.
         END IF
!
         energy2 = energy4(i,i)
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if(k > 0)then
               energy2 = energy2 + energy4(k,k)
            end if
         end do
         energy2 = energy2 + fscale*repul_energy(i,i)
!
         del_energy = energy2 - energy1

         if(special_atom)then
            if(quench)then
               EEE = del_energy/E0
            else
               EEE = del_energy/delt
            end if
         else
            EEE = del_energy/Etot
         end if

         IF(EEE <  0.0_wp) goto 22
         IF(EEE > 50.0_wp) goto 11
         IF(rand() < exp(-EEE)) goto 22
11       CONTINUE
            rxyz(i,1:3) = coprxyz(1:3)
            nfail = nfail + 1
            if(update)then
               call new_nlist
            end if
22       CONTINUE
!        IF (CELL(rxyz(i,:)) /= CELL(coprxyz)) call NEW_NLIST
      END DO

      END DO
!     print '(f0.6,2(i7,1x),f0.6)',de,NR*(ilast-ifirst+1),nfail,real(nfail)/(NR*(ilast-ifirst+1))
!     if(real(nfail)/(NR*(ilast-ifirst+1)) < 0.5)then
!        de = de*1.01
!     else
!        de = de/1.01
!     end if
   END SUBROUTINE RELAX

END MODULE  relax_mod

