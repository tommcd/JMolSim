
MODULE VERLET_LIST_MOD
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom_max,natom
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nabormx = 100
   real(wp),allocatable,private:: rxyz_old(:,:)
   real(wp),private:: rv,rv2,skin,skin2,rcut
   integer,allocatable:: nlist(:),list(:,:)
   integer,allocatable:: cnlist(:),clist(:,:)
   integer:: cclst(nabormx)
!
!  rxyz_old(i,3) : position of particle i at the time the list was made
!  rv              : Verlet radius
!  rv2=rv*rv
!  skin =rv-rc     : difference Verlet radius and cut-off
!  skin2=skin*skin
!  nlist(i)        : number of particles in Verlet list of particle i
!  list(i,j)       : Verlet-list of particle i
!
CONTAINS

   SUBROUTINE INIT_VERLET_LIST(rv0,rcut0)
      real(wp),intent(in):: rv0,rcut0
      rv = rv0
      rcut = rcut0
      RV2 = RV*RV
      SKIN = RV - rcut
      SKIN2 = SKIN*SKIN
      allocate( rxyz_old(natom_max,3) )
      allocate( nlist(natom_max),list(natom_max,nabormx) )
      allocate( cnlist(natom_max),clist(natom_max,nabormx) )
   END SUBROUTINE

   SUBROUTINE NEW_VLIST
!     makes the Verlet list using the neighbour list
      USE nlist_mod
!!USE sort_mod
      real(wp):: r(3),r2
      integer:: i,j,ic,jj,nc
      nlist(1:natom) = 0
      rxyz_old(1:natom,1:3) = rxyz(1:natom,1:3)
      call NEW_NLIST(1,natom)
      do i = 1, natom
      ic = CELL(rxyz(i,1:3))  ! link list
      call NEIGCELL(ic,1,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         j = hoc(nc)
         do while(j /= 0)
            if (j /= i .AND. j > i) then
            r(1:3) = rxyz(j,1:3)- rxyz(i,1:3)
            call pbc(r)
            r2 = r(1)**2 + r(2)**2 + r(3)**2
            if (r2 < RV2) then
               NLIST(i) = NLIST(i) + 1
               NLIST(j) = NLIST(j) + 1
               LIST(i, NLIST(i)) = j
               LIST(j, NLIST(j)) = i
            end if
            end if
            j = LL(j)
         end do
      end do cell_loop
      end do
      if (ANY(NLIST(1:natom) >= nabormx) ) stop 'NLIST(:) >= nabormx'
!!do i = 1,natom
!!   call qsort(nlist(i),list(i,1:nlist(i)),cclst)
!!   list(i,1:nlist(i)) =cclst(1:nlist(i))
!!end do
   END SUBROUTINE NEW_VLIST

   logical FUNCTION update(i)
      integer,intent(in):: i
      real(wp):: r(3)
      r(1:3) = rxyz(i,1:3) - rxyz_old(i,1:3)
      call pbc(r)
      update = (dot_product(r,r) > SKIN2*0.25_wp)
   END FUNCTION

   FUNCTION UPDATE2(ifirst,ilast)
!     decides whether the list needs to be reconstructed
      logical:: UPDATE2
      integer,intent(in):: ifirst,ilast
      real(wp):: dispmx
!     a conservative test of the list skin crossing
      dispmx = MAXVAL(ABS( rxyz(ifirst:ilast,1:3) - rxyz_old(ifirst:ilast,1:3) ))
      dispmx = 2.0_wp*SQRT( 3.0_wp*dispmx**2 )
      UPDATE2 = ( dispmx > skin )
   END FUNCTION

END MODULE VERLET_LIST_MOD

