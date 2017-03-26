
MODULE Keating
   USE precision_mod, only: wp
   USE Keating_parameters
   implicit none
   PUBLIC energy4
CONTAINS

   PURE FUNCTION energy4(ifirst,ilast)
!-----function for the Keating energy
      USE coordinates_mod
      USE atom_types, only: atom,ncmax,iHydrogen
      USE connectivity_mod, only: proximity

      USE matvec3d_mod,only: dotp_3d,len_3d
      real(wp):: energy4
      integer,intent(in):: ifirst,ilast
      real(wp):: ebond,en1,eangle,rij(3),rik(3) !,csa
      integer:: i,L,M,j,k
      real(wp):: kang(0:2,0:2,0:2)
      kang(0,1,0)= KOSiO
      kang(2,1,0)= KOSiO
      kang(0,1,2)= KOSiO
      kang(2,1,2)= KOSiO
      !! kang(0,1,1)= KSiSiO  ! Not used for now (no Si-Si bonds)
      !! kang(1,1,0)= KSiSiO
      !! kang(1,1,1)= KSiSiSi
      kang(1,0,1)= KSiOSi
      kang(1,2,1)= KSiOSi
      ebond = 0.0_wp
      eangle = 0.0_wp

!
print*, 'ebond =',ebond
      do i = ifirst,ilast

print*, 'i (keating) = ', i
      if (atom(i) == iHydrogen) cycle
!     bond energy
      do L = 1,ncmax(atom(i))
         j = proximity(i,L)
         if (atom(j) == iHydrogen) cycle
         if (j /= 0) then
         rij(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
         call pbc(rij)
         ebond = ebond + 0.5_wp*kbond(atom(i)+atom(j)) &
               *(sqrt(dot_product(rij,rij)) - abond(atom(i)+atom(j)))**2
!!             *(len_3d(rij) - abond(atom(i)+atom(j)))**2
         end if
      end do
!
!     angle energy
      do L = 1,ncmax(atom(i))-1
         j = proximity(i,L)
         if (atom(j) == iHydrogen) cycle
         if (j /= 0) then
         do M = L + 1,ncmax(atom(i))
            k = proximity(i,M)
            if (atom(k) == iHydrogen) cycle
            if (k /= 0) then
            rij(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
            rik(1:3) = rxyz(k,1:3) - rxyz(i,1:3)
            call pbc(rij)
            call pbc(rik)
!!          csa = dot_product(rij,rik)/(sqrt(dot_product(rij,rij))*sqrt(dot_product(rik,rik)))
!!          en1 = (csa + acosang(atom(i)))**2
            en1 = (dotp_3d(rij,rik)/(len_3d(rij)*len_3d(rik)) + acosang(atom(i)))**2
            eangle = eangle + 0.5_wp*kang(atom(j),atom(i),atom(k))*en1
         end if
         end do
         end if
      end do

      end do
      energy4 = ebond + eangle
      RETURN
    END FUNCTION energy4

END MODULE Keating

