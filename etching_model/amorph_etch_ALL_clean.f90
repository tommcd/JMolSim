!!>include 'precision_mod.f90'
MODULE precision_mod
   implicit none
!  integer, parameter :: sp = kind(1.0)
!  integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: sp = selected_real_kind(6,30)
   integer, parameter :: dp = selected_real_kind(15,300)
   integer, parameter :: qp_preferred = selected_real_kind(30,1000)
   integer, parameter :: qp = (1+sign(1,qp_preferred))/2*qp_preferred+ &
                              (1-sign(1,qp_preferred))/2*dp
!
   integer,parameter,public :: wp = dp
   integer,parameter,public :: i4b = selected_int_kind(9)
END MODULE precision_mod

!!>include 'global_vars_mod.f90'
MODULE global_vars_mod
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 1.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, etot
   integer:: natom,nseed,natom_max
   integer:: nattached, nrelax
END MODULE global_vars_mod

!!>include 'files_mod.f90'
MODULE files_mod
   public get_free_file_unit, myflush
   contains

   SUBROUTINE get_free_file_unit(iu)
      implicit none
      character(len=*),parameter :: sub_name='get_free_file_unit'
      integer,intent(out):: iu
      integer,parameter:: istart=11
      integer,parameter:: unit_max=100000
      integer:: ios
      logical:: lopen,lexists
      do iu = istart,unit_max
         inquire(unit=iu, iostat=ios, exist=lexists, opened=lopen)
         if (ios == 0) then
            if (lexists .and. .not.lopen) return
         else
            write(*,*) sub_name,': error, iostat = ',ios
            stop
         end if
      end do
      write(*,*) sub_name,': error, no free units'
      stop
   END SUBROUTINE get_free_file_unit

   SUBROUTINE myflush(iu)
      integer,intent(in):: iu
      character(132):: fn
      inquire(unit=iu,name=fn)
      close(iu)
      open(unit=iu,file=trim(fn),position='append')
   END SUBROUTINE myflush

END MODULE files_mod

!!>include 'sort_mod.f90'

MODULE sort_mod
      implicit none
! Quicksort modified from Numerical Recipes
CONTAINS

   SUBROUTINE qsort(n,arr0,arr)
      integer,intent(in):: n,arr0(:)
      integer,intent(out):: arr(:)
      integer,parameter:: M=7,NSTACK=50
      integer:: i,ir,j,jstack,k,L,istack(NSTACK)
      integer:: a,temp
      jstack=0
      L=1
      ir=n
      arr(1:n)=arr0(1:n)
1     if (ir-L < M) then
        do j=L+1,ir
          a=arr(j)
          do i=j-1,1,-1
            if (arr(i) <= a) GOTO 2
            arr(i+1)=arr(i)
          end do
          i=0
2         arr(i+1)=a
        end do
        if (jstack == 0)return
        ir=istack(jstack)
        L=istack(jstack-1)
        jstack=jstack-2
      else
        k=(L+ir)/2
        temp=arr(k)
        arr(k)=arr(L+1)
        arr(L+1)=temp
        if (arr(L+1) > arr(ir)) then
          temp=arr(L+1)
          arr(L+1)=arr(ir)
          arr(ir)=temp
        end if
        if (arr(L) > arr(ir)) then
          temp=arr(L)
          arr(L)=arr(ir)
          arr(ir)=temp
        end if
        if (arr(L+1) > arr(L)) then
          temp=arr(L+1)
          arr(L+1)=arr(L)
          arr(L)=temp
        end if
        i=L+1
        j=ir
        a=arr(L)
3       continue
          i=i+1
        if (arr(i) < a) GOTO 3
4       continue
          j=j-1
        if (arr(j) > a) GOTO 4
        if (j < i) GOTO 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        GOTO 3
5       arr(L)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if (jstack > NSTACK) stop 'NSTACK too small in qsort'
        if (ir-i+1 >= j-L) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=L
          L=i
        end if
      end if
      GOTO 1
   END SUBROUTINE
!  (C) Copr. 1986-92 Numerical Recipes Software #1-0zV'n26)B3.


   SUBROUTINE shell(n,v)
      integer,intent(in):: n
      integer,intent(inout):: v(:)
! Sorts vector v(1:n) into ascending numerical order
! by Shell's method (diminishing increment sort)
      integer:: i,j,inc,b
      inc = 1   ! Determine the starting increment
1     inc = 3*inc+1
      if (inc <= n) GOTO 1
2     continue   ! Loop over the partial sorts
      inc = inc/3
      do i=inc+1,n   ! Outer loop of straight insertion.
         b=v(i)
         j=i
3        if (v(j-inc) > b) then   ! Inner loop of straight insertion.
            v(j)=v(j-inc)
            j=j-inc
            if (j <= inc) GOTO 4
            GOTO 3
         end if
4        v(j)=b
      end do
      if (inc > 1) GOTO 2
      return
   END SUBROUTINE

   SUBROUTINE sort3(iv)
      integer,intent(inout):: iv(:)
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
      if (iv(3) < iv(2)) call swap(iv(3),iv(2))
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
   contains
      SUBROUTINE swap(x, y)
         integer,intent(inout):: x,y
         integer:: tmp
         tmp = x
         x = y
         y = tmp
      END SUBROUTINE
   END SUBROUTINE

END MODULE

!!>include 'atom_types_mod.f90'

MODULE atom_types_mod
   USE precision_mod, only: wp
   implicit none
   integer,parameter:: iOxygen = 0
   integer,parameter:: iSilicon = 1
   integer,parameter:: iOxygenH = 2
   integer,parameter:: iHydrogen = 3
   integer,parameter:: iOw=4
   integer,parameter:: iHw=5
   integer,parameter:: ntyp = 5
   character(2),parameter:: atom_name(0:ntyp)=(/ 'O ','Si','OH','H ','Ow','Hw' /)
   integer,parameter:: ncmax(0:ntyp)=(/2,4,2,1,2,1/)
   integer,allocatable:: atom(:)
CONTAINS

   pure function name2atom(c)
      integer:: name2atom
      character(*),intent(in):: c
      integer:: i
      do i = 0,ntyp
         if(trim(c) == trim(atom_name(i)))then
            name2atom = i
            exit
         end if
      end do
      if (i > ntyp) name2atom = -1
   end function

END MODULE atom_types_mod

!!>include 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:)
   real(wp):: boxl(3),boxli(3),boxl2(3) !,boxl2i,boxl2n
   interface pbc
      module procedure  pbc_v, pbc_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
if(abs(r(1)) > boxl2(1)) r(1) = r(1) - sign(boxl(1),r(1))
if(abs(r(2)) > boxl2(2)) r(2) = r(2) - sign(boxl(2),r(2))
if(abs(r(3)) > boxl2(3)) r(3) = r(3) - sign(boxl(3),r(3))
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      integer:: i
      do i=1,size(r,1)
         if(abs(r(i,1)) > boxl2(1)) r(i,1) = r(i,1) - sign(boxl(1),r(i,1))
         if(abs(r(i,2)) > boxl2(2)) r(i,2) = r(i,2) - sign(boxl(2),r(i,2))
         if(abs(r(i,3)) > boxl2(3)) r(i,3) = r(i,3) - sign(boxl(3),r(i,3))
      end do
   END SUBROUTINE

   SUBROUTINE WRITE_XYZ(iu,nattached)
      USE atom_types_mod, only: atom_name,atom
      integer,intent(in):: iu,nattached
      integer:: i
      write(iu,*) natom
      write(iu,'(a,i6)') 'Silica   nattached = ',nattached
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*Angstrom
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,it
      real(wp):: x,y,z
      character(32):: ctmp
      read(iu,*) natom
      read(iu,*) ctmp!,ctmp,ctmp,it
      it = 0
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod

!!>include 'connectivity_mod.f90'

MODULE connectivity_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,nseed
   USE atom_types_mod
   implicit none
   integer,allocatable:: proximity(:,:)
!
CONTAINS

   SUBROUTINE write_proximity(iu,nattached)
!     write topology/connectivity information to file unit iu
      integer,intent(in):: iu,nattached
      integer:: i,j,k
      write(iu,*)'# nattached = ',nattached
      write(iu,*) natom
      do i = 1,natom
         k = count(proximity(i,:) > 0)
         write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
            i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE write_proximity


   SUBROUTINE read_proximity(iu,nattached)
!     Read topology/connectivity information from file unit iu
      integer,intent(in):: iu,nattached
      integer:: i,j,ic,ii,it
      character(32):: ctmp
      read(iu,*) ctmp,ctmp,ctmp,it
      !if (it /= nattached)stop 'it /= nattached'
      read(iu,*) natom
      do i = 1,natom
         read(iu,*) ii,ctmp,ic,(proximity(i,j),j = 1,ncmax(name2atom(ctmp)))
         atom(i) = name2atom(ctmp)
         if ( atom(i) < 0 ) stop 'error reading atom type'
      end do
      rewind(iu)
   END SUBROUTINE read_proximity


   SUBROUTINE set_proximity(iat,iold,inew)
!     Where the value of proximity is iold set it to inew.
!     Used to break/switch bonds.
!     e.g. the following move requires two call to set_proximity
!
!     A--B         A  B      set_proximity(A,B,C)  switch A-B to A-C
!            -->    \        set_proximity(B,A,0)  break B-A
!        C           C
!
      integer,intent(in):: iat,iold,inew
      integer:: i
      do i = 1,ncmax(atom(iat))
         if (proximity(iat,i) == iold) then
            proximity(iat,i) = inew
            exit
         end if
         if (i == ncmax(atom(iat))) STOP 'set_proximity: ERROR'
      end do
   END SUBROUTINE set_proximity


   PURE FUNCTION nearest_neighbor2(L,M)
!     Are atoms L & M nearest neighbors up to order 2?
!     i.e. Are they bonded together or to a common atom K
!                   K
!     L--M   or    / \
!                 L   M
      logical:: nearest_neighbor2
      integer,intent(in):: L,M
      integer:: i,j,K
      nearest_neighbor2 = .false.
      do i = 1,ncmax(atom(L))
         K = proximity(L,i)
         if (K == 0) cycle
         if (K == M) then
            nearest_neighbor2 = .true.
            return
         end if
         do j = 1,ncmax(atom(K))
            if (proximity(K,j) == M) then
               nearest_neighbor2 = .true.
               return
            end if
         end do
      end do
   END FUNCTION nearest_neighbor2


   PURE FUNCTION is_dangling_group(iSi,iat)
!     Is the Silicon (iSi) attached to bridging Oxygen (iat)
!     in a 'dangling' group? i.e.  --O  <--iat
!                                     \
!     is it bonded to 3(OH)       HO--Si--OH
!     groups?                          |
!                                      OH
!
      integer,intent(in):: iSi,iat
      logical:: is_dangling_group
      integer:: k,j
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == iat) cycle
         if (j == 0) cycle
         if (atom(j) /= iOxygenH) then
            is_dangling_group = .FALSE.
            RETURN
         end if
      end do
      is_dangling_group = .TRUE.
   END FUNCTION


   PURE FUNCTION noh(i)
! number of OH atoms bonded to atom i
! is count(atom(proximity(i,1:ncmax(atom(i)))) == iOxygenH)
      integer,intent(in):: i
      integer:: noh,j,k
      noh = 0
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         if(atom(k) == iOxygenH) noh = noh + 1
      end do
   END FUNCTION


   SUBROUTINE check_proximity(consistent,ifirst,k)
!     Check that if atom i is bonded to k then
!     atom k is bonded to i, if everything is ok
!     then consistent is .true.
      logical,intent(out):: consistent
      integer,intent(out):: ifirst,k
      integer:: i,j
      consistent = .true.
      ifirst = 0
      do i = 1,natom
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (k == 0) cycle
            if ( count(proximity(k,:) == i) /= 1) then
               consistent = .false.
               ifirst = i
               return
            end if
         end do
      end do
   END SUBROUTINE

   SUBROUTINE check_proximity2(consistent,ifirst,k)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      logical,intent(out):: consistent
      integer,intent(out):: ifirst,k
      integer:: i,j,ia,a1,a2,i1,i2
      consistent = .true.
      ifirst = 0
      do i = 1,natom
         ia = atom(i)
         select case(ia)
         case(iHydrogen)  ! Check H is bonded to 1 OH Oxygen
            k = proximity(i,1)
            if (k == 0) GOTO 100
            if (atom(k) /= iOxygenH) GOTO 100
            if (any(proximity(i,2:) /= 0)) GOTO 100
         case(iSilicon)  ! Check Si is bonded to 4 atoms
            do j = 1,ncmax(ia)
               k = proximity(i,j)
               if (k == 0) GOTO 100
               if (atom(k)==iSilicon) GOTO 100   ! but not Si
               if (atom(k)==iHydrogen) GOTO 100  ! and not H
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = proximity(i,1)
            k = proximity(i,1)
            if (i1 == 0) GOTO 100
            a1 = atom(proximity(i,1))
            if (a1 /= iSilicon) GOTO 100
            if (i > nseed) then  ! non seed bridging O
               i2 = proximity(i,2)  !  must be bonded to  2 Silicons
               k = proximity(i,2)
               if (i2 == 0) GOTO 100
               a2 = atom(proximity(i,2))
               if (a2 /= iSilicon) GOTO 100
            end if
            if (any(proximity(i,3:) /= 0)) GOTO 100
!         case(iOxygenH)  ! OH Oxygen must be bonded to Si and H
!            i1 = proximity(i,1)
!            k = proximity(i,1)
!            if (i1 == 0) GOTO 100
!            a1 = atom(proximity(i,1))
!            if (i > nseed) then
!               i2 = proximity(i,2)
!               k = proximity(i,2)
!               if (i2 == 0) GOTO 100
!               a2 = atom(proximity(i,2))
!               if (i2 == 0) GOTO 100
!               if (.NOT.((a1 == iHydrogen .and. a2 == iSilicon).or. &
!                         (a2 == iHydrogen .and. a1 == iSilicon))) GOTO 100
!            else
!               if (a1 /= iHydrogen) GOTO 100
!            end if
!            if (any(proximity(i,3:) /= 0)) GOTO 100

         case(iOxygenH)  ! OH combined oxygen-hydrogen
            k = proximity(i,1)
            if (k==0) GOTO 100
            if (atom(k)/=iSilicon) GOTO 100
            if (any(proximity(i,2:)/=0)) GOTO 100
         case default
            stop 'unknown atom type'
         end select
      end do
      return
100   consistent = .false.
      ifirst = i
      return
   END SUBROUTINE

   SUBROUTINE delete_atom(i)
!     Delete atom i from the system and renumber the
!     last atom (natom) to number i
      USE coordinates_mod, only: rxyz
      integer,intent(in):: i
      integer:: j,k
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0)cycle
         where( proximity(k,:) == i ) proximity(k,:) = 0
      end do
      rxyz(i,1:3) = rxyz(natom,1:3)
      proximity(i,:) = proximity(natom,:)
      atom(i) = atom(natom)
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0)cycle
         !print *,'i = ',i
         !print *,'j = ',j
         !print *,'k = ',k
         where( proximity(k,:) == natom ) proximity(k,:) = i
      end do
      rxyz(natom,1:3) = 0.0_wp
      proximity(natom,:) = 0
      atom(natom) = 0
      natom = natom - 1
   END SUBROUTINE delete_atom

   SUBROUTINE delete_group0(iSi,iat)
!     delete 'dangling' SiO3 group attached to
!     bridging Oxygen iat
      USE sort_mod,only:shell
      integer,intent(in):: iSi,iat
      integer:: i,k,j,ilst(4),ni
      ni = 0
      ilst = 0
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == iat) cycle
         if (j == 0) cycle
         if (atom(j) /= iOxygenH) then
            write(*,*)'delete_group: error'
            stop 'not a dangling group'
         end if
         ni = ni+1
         ilst(ni) = j
      end do
      ni = ni+1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE

   SUBROUTINE delete_group(iSi)
!     Delete Silicon iSi and all 4 Oxygens bonded to it
!     Delete atom iSi and its nearest neighbors
      USE sort_mod,only:shell
      integer,intent(in):: iSi
      integer:: i,k,j,ilst(9),ni
      ni = 0
      ilst = 0
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == 0) cycle
         if (atom(j) /= iOxygenH) then
            write(*,*)'delete_group: error'
            stop 'not a dangling group'
         end if
         ni = ni+1
         ilst(ni) = j
      end do
      ni = ni+1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE

   SUBROUTINE delete_group2(iSi)
!     Delete Si(OH)4 cluster
!     Delete atom iSi and its nearest neighbors up to order 2
      USE sort_mod,only:shell
      integer,intent(in):: iSi
      integer:: k,j,i,m
      integer:: ilst(9),ni
      ni = 0
      ilst = 0
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == 0) cycle
         if (atom(j) /= iOxygenH) then
            write(*,*)'delete_group: error'
            stop 'not a dangling group'
         end if
         do m = 1,ncmax(atom(j))
            i = proximity(j,m)
            if (i == iSi) cycle
            if (i == 0) cycle
            if (atom(i) /= iHydrogen) then  !only for debugging
               STOP 'ERROR: delete_group2'
            end if
            ni = ni+1
            ilst(ni) = i
         end do
         ni = ni+1
         ilst(ni) = j
      end do
      ni = ni+1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE

END MODULE connectivity_mod

!!>include 'nlist_mod_2.f90'
MODULE NLIST_MOD
      USE precision_mod, only: wp
      USE global_vars_mod, only: natom_max,natom
      USE coordinates_mod, only: boxl,rxyz,boxl2
      implicit none
      PRIVATE ! assigns every variable as private unless otherwise stated
      integer,public:: neigh
      real(wp),private:: delx,dely,delz,delmin
      real(wp),private:: delxi,delyi,delzi
      integer:: ncelx,ncely,ncelz,ncelt
      integer,parameter:: ncelmax=100000,neighmx=150
      integer,allocatable,public:: hoc(:),hoc_old(:)
      integer,public::ncell(neighmx)
      integer,allocatable,public:: ll_old(:),ll(:),lr(:)
      ! when use the global PRIVATE it is required to specify whether each of
      ! the subroutines is of public type if that is desired
      public:: NEIGCELL,INIT_NLIST,RE_INIT_NLIST,NEW_NLIST,print_cell!,CELL_OB,CELL_OC
      public:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll,cell
!
!     ll(i)      : linked list particle i
!     hoc(ic)    : head of chain cell ic
!     ncelx,y,z  : number of cells in x, y or z direction
!     ncelt      : total number of cells
!     ncelmax    : maximum number of cells (change and
!                  recompile if larger number is needed)
!     neigh      : number of cells for interactions
!
CONTAINS

! the following function is not the first thing called so it can be put anywhere
! in the file really. INIT_NLIST() is the first section called
   pure function nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   end function

   SUBROUTINE INIT_NLIST(boxl,Rc)
      real(wp),intent(in):: Rc,boxl(:)
      real(wp)::Lx,Ly,Lz
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max))
      Lx = boxl(1)
      Ly = boxl(2)
      Lz = boxl(3)
! ncel_ variables give no. of cells in respective directions
      ncelx = INT(Lx/Rc)
      ncely = INT(Ly/Rc)
      ncelz = INT(Lz/Rc)
      ncelt = ncelx*ncely*ncelz
      allocate(hoc(ncelt),hoc_old(ncelt))
! del_ variables specify the actual length of each cell in respective direction
      delx = Lx/ncelx
      dely = Ly/ncely
      delz = Lz/ncelz
! delmin is the smallest length of a cell side comparing x,y and z directions
      delmin = min(delx,dely,delz)
! del_i is the inverse of repestive cell lengths
      delxi = 1.0_wp/delx
      delyi = 1.0_wp/dely
      delzi = 1.0_wp/delz
!      write(*,*)'ncelx = ',ncelx
!      write(*,*)'ncely = ',ncely
!      write(*,*)'ncelz = ',ncelz
!      write(*,*)'ncelt = ',ncelt
!      write(*,*)'delx = ',delx
!      write(*,*)'dely = ',dely
!      write(*,*)'delz = ',delz
      LL = 0
      HOC = 0
      LR = 0 ! although it is not really used, LR is essentially a reverse link
             ! list used to go back through elements
   END SUBROUTINE INIT_NLIST

   SUBROUTINE RE_INIT_NLIST(boxl,Rc)
      real(wp),intent(in):: Rc,boxl(:)
      real(wp)::Lx,Ly,Lz
      if (allocated(ll_old)) then
         deallocate(ll_old,ll,lr)
         deallocate(hoc,hoc_old)
      end if
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max))
      Lx = boxl(1)
      Ly = boxl(2)
      Lz = boxl(3)
      ncelx = INT(Lx/Rc)
      ncely = INT(Ly/Rc)
      ncelz = INT(Lz/Rc)
      ncelt = ncelx*ncely*ncelz
      ncelx = INT(Lx/Rc)
      ncely = INT(Ly/Rc)
      ncelz = INT(Lz/Rc)
      ncelt = ncelx*ncely*ncelz
      allocate(hoc(ncelt),hoc_old(ncelt))
      delx = Lx/ncelx
      dely = Ly/ncely
      delz = Lz/ncelz
      delmin = min(delx,dely,delz)
      delxi = 1.0_wp/delx
      delyi = 1.0_wp/dely
      delzi = 1.0_wp/delz
      LL = 0
      HOC = 0
      LR = 0
   END SUBROUTINE RE_INIT_NLIST

   PURE FUNCTION CELL(XYZ)
!     determines cell number for position x,y,z
!     with origin at centre of box
      real(wp),intent(in)::  XYZ(:)
      integer:: CELL
      CELL = INT((XYZ(1)+boxl2(1))*delxi) &
           + INT((XYZ(2)+boxl2(2))*delyi)*ncelx &
           + INT((XYZ(3)+boxl2(3))*delzi)*ncelx*ncely + 1
!print *,'xyz ',xyz
!print *,boxl2
!print *,delxi,delyi,delzi
!print *,'cell ',cell
!if(cell < 0)stop
      RETURN
   END FUNCTION CELL

!   PURE FUNCTION CELL_OC(XYZ)
!!     determines cell number for position x,y,z
!!     with origin  at centre of box
!      real(wp),intent(in)::  XYZ(:)
!      integer:: CELL_OC
!      CELL_OC = INT((XYZ(1)+boxl2)*delxi) &
!              + INT((XYZ(2)+boxl2)*delyi)*ncelx &
!              + INT((XYZ(3)+boxl2)*delzi)*ncelx*ncely + 1
!      RETURN
!   END FUNCTION CELL_OC

!   PURE FUNCTION CELL_OB(XYZ)
!!     determines cell number for position x,y,z
!!     with origin at bottom of box
!      real(wp),intent(in)::  XYZ(:)
!      integer:: CELL_OB
!      CELL_OB = INT((XYZ(1)+boxl2)*delxi) &
!              + INT((XYZ(2)+boxl2)*delyi)*ncelx &
!              + INT( XYZ(3)       *delzi)*ncelx*ncely + 1
!      RETURN
!   END FUNCTION CELL_OB


   SUBROUTINE NEW_NLIST
!     makes a new neighbour list using the linked-list algorithm
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head-of-chain
      ! make linked list:
      do i = 1,natom
         ! determine cell number
         ic = CELL(rxyz(i,:))
!print*,'i = ',i
if (ic>ncelt) then
print*,'rxyz(',i,',:) = ',rxyz(i,:)
print*,'ic = ',ic
end if
!print*,'LL(i)= ',LL(i)
!print*,'HOC(ic) = ',HOC(ic)
!print*,''
         ! update linked-list and head of chain
         LL(i) = HOC(ic)
         if (HOC(ic)/=0)LR(HOC(ic)) = i
         HOC(ic) = i
      end do
      RETURN
   END SUBROUTINE NEW_NLIST

   SUBROUTINE PUSH_CELL(i,ic)
      integer,intent(in):: i,ic
      LL(i) = HOC(ic)
      if (HOC(ic)/=0)LR(HOC(ic)) = i
      HOC(ic) = i
   END SUBROUTINE

   SUBROUTINE POP_CELL(i,ic)
      integer,intent(in):: i,ic
      if (HOC(ic) == i) then
         HOC(ic) = LL(i)
         if (LL(i) /= 0)LR(LL(i)) = 0
         LL(i) = 0
      else
         LL(LR(i)) = LL(i)
         if (LL(i) /= 0)LR(LL(i)) = LR(i)
         LL(i) = 0
      end if
   END SUBROUTINE

   PURE SUBROUTINE NEIGCELL(ic,nlayer,neigh,ncell)
!     determines the neigh neighbours
! NOTE: NOT Periodic in Z-direction
      integer,intent(in):: ic,nlayer
      integer,intent(out):: neigh
      integer,intent(out):: ncell(:)
      integer ix,iy,iz,inn,icx,icy,icz,iccx,iccy,iccz,n2
      neigh = 0
      n2 = ncelx*ncely
      icz = ic/n2 + 1
      if (mod(ic,n2)==0)icz = icz - 1
      icy = (ic-(icz-1)*n2)/ncelx + 1
      if (mod(ic-(icz-1)*n2,ncelx)==0)icy = icy - 1
      icx = ic-(icy-1)*ncelx-(icz-1)*n2
      ncell = 0
      do iz = -nlayer,nlayer
         iccz = icz + iz
         if (iccz < 1) then
            iccz = iccz + ncelz ! Periodic in Z-direction
         else if (iccz > ncelz) then
            iccz = iccz - ncelz  ! Periodic in Z-direction
         end if
         do iy = -nlayer,nlayer
            iccy = icy + iy
            if (iccy < 1) then
               iccy = iccy + ncely
            else if (iccy > ncely) then
               iccy = iccy - ncely
            end if
            do ix = -nlayer,nlayer
               iccx = icx + ix
               if (iccx < 1) then
                  iccx = iccx + ncelx
               else if (iccx > ncelx) then
                  iccx = iccx - ncelx
               end if
               inn = iccx + (iccy-1)*ncelx + (iccz-1)*n2
               neigh = neigh + 1
               ncell(neigh) = inn
            end do
         end do
      end do
      RETURN
   END SUBROUTINE NEIGCELL

!   SUBROUTINE NEIGCELL_B(ic)
      !! not used (unfinished)
      !integer,intent(in):: ic
      !integer:: icx,icy,icz,n2
      !n2 = ncelx*ncely
      !icz = ic/n2 + 1
      !icy = (ic-(icz-1)*n2)/ncelx + 1
      !icx = ic-(icy-1)*ncelx-(icz-1)*n2
      !neigh = 27
      !ncell(13) = ic - 1
      !ncell(14) = ic
      !ncell(15) = ic + 1
      !ncell(10:12) = ncell(13:15) - ncelx
      !ncell(16:18) = ncell(13:15) + ncelx
      !ncell( 1: 9) = ncell(10:18) - n2
      !ncell(19:27) = ncell(10:18) + n2
   !END SUBROUTINE

   SUBROUTINE print_cell(ic,iu)
      integer,intent(in):: ic,iu
      integer:: icx,icy,icz
      icz = ic/(ncelx*ncely) + 1
      if (mod(ic,(ncelx*ncely))==0)icz = icz - 1
      icy = (ic-(icz-1)*(ncelx*ncely))/ncelx + 1
      if (mod(ic-(icz-1)*(ncelx*ncely),ncelx)==0)icy = icy - 1
      icx = ic-(icy-1)*ncelx-(icz-1)*(ncelx*ncely)
      write(iu,*) 'ic = ',ic
      write(iu,*) 'icx = ',icx
      write(iu,*) 'icy = ',icy
      write(iu,*) 'icz = ',icz
   END SUBROUTINE

!EXAMPLE OF USE OF LINK LIST
!!     ---determine cell number
!      ic = CELL(rxyz(i,:))
!!     ---determine neighbour cells
!      call NEIGCELL(ic)
!!        ---loop over neighbours and same cell
!         do inn = 1, neigh
!            nc = ncell(inn)
!            if (nc == 0)cycle
!            j = HOC(nc)
!            do while (j /= 0)
!               if (j /= I) then
!                  ! Calc energy etc.
!               end if
!!              ---next particle in the cell
!               j = LL(j)
!            end do
!         end do

   PURE SUBROUTINE Z_NEIGH_CELL(iz,neigh,ncell)
!     determines the neigh cells in the iz'th z-layer
      integer,intent(in):: iz
      integer,intent(out):: neigh
      integer,intent(out):: ncell(:)
      integer:: n2,i
      n2 = ncelx*ncely
      neigh = n2
      forall(i=1:n2) ncell(i) = n2*(iz-1) + i
   END SUBROUTINE


END MODULE NLIST_MOD

!!>include 'HKNonLattice2_mod.f90'
MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster=0,n_cluster_old=0
   integer,parameter,private :: ncmx = 2000
   integer,allocatable:: atomL(:),atomL_old(:)
!   integer,allocatable,private:: cluster(:,:)
   integer,private:: clusterC(ncmx)
CONTAINS

   SUBROUTINE Init_HKNonLattice(natom_max)
      integer,intent(in):: natom_max
      if(allocated(atomL))then
         deallocate(atomL,atomL_old)
         allocate(atomL(natom_max),atomL_old(natom_max))
!         allocate(cluster(ncmx,natom_max))
      else
         allocate(atomL(natom_max),atomL_old(natom_max))
      end if
      atomL=0
      atomL_old=0
   END SUBROUTINE Init_HKNonLattice

   SUBROUTINE HKNonLattice(NumberOfNodes,NodeNext,NumberOfClusters,NodeL)
      USE sort_mod, only: qsort
!===============================================================
!     Adaptation of Hoshen--Kopelman cluster labeling algorithm
!     for molecular networks.
!     A simplified version of the Algorithm by
!     Ahmed AL-Futaisi and Tadeusz Patzek
!     Translated to Fortran 95: T.C.McDermott 14/06/05
!===============================================================
! Input arguments:
!     NumberOfNodes = Number of nodes (atoms) in network (molecule)
!     NodeNext = Neighboring nodes connected to each node
      integer,intent(in):: NumberOfNodes, NodeNext(:,:)
!
! Output arguments:
!     NumberOfClusters = Number of occupied clusters
!     NodeL = Cluster labels of nodes
      integer,intent(out):: NumberOfClusters, NodeL(:)
!
      integer:: NodeLP(NumberOfNodes),NodeLP1(NumberOfNodes),iNodeLP,iCluster
      integer:: N(size(NodeNext,dim=2)),NodeNextL(size(NodeNext,dim=2))
      integer:: i,ii,j,nnmax,iNodeNextL,inn,NodeLPmin,k,m
      integer:: RelabL1(NumberOfNodes),RelabL(NumberOfNodes),RelabLB(NumberOfNodes)
      !
      ! STEPS 1,2 & 3 of AL-Futaisi and Tadeusz Patzek algorithm
      !
      nnmax=size(NodeNext,dim=2)
      NumberOfClusters = 0
      NodeL = 0
      iNodeLP = 0
      NodeLP = 0   ! Array used for relabeling steps
      iCluster = 0  ! iCluster counter
      !
      ! STEP 4: SCAN THE NETWORK NODES
      !
      do i = 1,NumberOfNodes
         if ( ANY(NodeNext(i,:) > 0) ) then
            iNodeNextL = 0
            NodeNextL = 0
            do ii = 1,nnmax
               j = NodeNext(i,ii)
               if (j == 0) cycle
               iNodeNextL = iNodeNextL + 1
               NodeNextL(iNodeNextL)=NodeL(j)
            end do
            if ( any(NodeNextL > 0) ) then  ! Case 4c ii: a labeled neighbor exists
               ! Put in the minimum labeling
               inn = 0
               do ii = 1,iNodeNextL
                  j = NodeNextL(ii)
                  if (j == 0) cycle
                  inn = inn + 1
                  N(inn) = j
               end do
               do k = 1,inn
                  M = NodeLP(N(k))
                  do while(M < N(k))
                     N(k) = M
                     M = NodeLP(N(k))
                  end do
               end do
               NodeLPmin = minval(N(1:inn))
               NodeL(i) = NodeLPmin
               NodeLP(N(1:inn)) = NodeLPmin
            else  ! Case 4c i: No labeled neighbour
               iCluster = iCluster + 1  ! Start a new cluster
               NodeL(i) = iCluster
               iNodeLP = iNodeLP + 1
               NodeLP(iNodeLP) = iCluster
            end if
         else ! This node is type 4b
            iCluster = iCluster + 1  ! Start a new cluster
            NodeL(i) = iCluster
            iNodeLP = iNodeLP + 1
            NodeLP(iNodeLP) = iCluster
         end if
      end do
      !
      ! STEP 5A: CORRECT LABELS IN NodeLP RECURSIVELY
      !
      do i = 1,iNodeLP
         k = i
         do while (NodeLP(k) < k)
            k = NodeLP(k)
         end do
         NodeLP(i) = k
      end do
      !
      ! STEP 5B: RENUMBER LABELS IN NodeLP TO RUN SEQUENTIALLY
      !
      call qsort(iNodeLP,NodeLP(1:iNodeLP),NodeLP1(1:iNodeLP))
      RelabLB(1:iNodeLP) = 0
      where(NodeLP1(2:iNodeLP) > NodeLP1(1:iNodeLP-1)) RelabLB(2:iNodeLP)=1
      RelabL1(1:iNodeLP-1) = NodeLP1(2:iNodeLP)*RelabLB(2:iNodeLP)
      RelabL = 0
      RelabL(1) = NodeLP1(1)
      ii = 1
      do i = 1,iNodeLP-1
         if (RelabL1(i)==0) cycle
         ii = ii + 1
         RelabL(ii) = RelabL1(i)
      end do
      do i = 1,ii
         where(NodeLP == RelabL(i)) NodeLP = i
      end do
      !
      ! STEP 6: APPLY THE CORRECT LABELS TO NodeL
      !
      do i = 1,iNodeLP
         where(NodeL == i) NodeL = NodeLP(i)
      end do
      !
      ! RECORD NUMBER OF CLUSTERS
      NumberOfClusters = maxval(NodeL)
      return
   END SUBROUTINE HKNonLattice




   SUBROUTINE cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
      integer,intent(in):: natom,n_cluster
      integer,intent(out):: maxcluster
      integer,allocatable,intent(out):: ClusterCount(:)
      integer:: i,L
      allocate(ClusterCount(n_cluster))
      ClusterCount(1:n_cluster) = 0
      do i = 1,natom
         L = AtomL(i)
         ClusterCount(L) = ClusterCount(L) + 1
      end do
      MaxCluster = maxloc(ClusterCount,1)
      ! some consistency checks
      do i = 1,n_cluster
         if (ClusterCount(i) /= count(atomL(1:natom) == i)) then
            stop 'error in analyse_cluster'
         end if
      end do
      if (sum(ClusterCount(1:n_cluster)) /= natom) then
         stop 'error in analyse_cluster: wrong count'
      end if
   END SUBROUTINE cluster_count






! some routines used mainly for debugging
!   SUBROUTINE analyse_cluster(natom)
!      implicit none
!      integer,intent(in):: natom
!      integer:: i,L
!      clusterC(1:n_cluster)=0
!!      cluster(1:n_cluster,1:ncmx)=0
!      cluster=0
!      do i = 1,natom
!         L = AtomL(i)
!         clusterC(L) = clusterC(L) + 1
!         cluster(L,clusterC(L))=i
!      end do
!      ! some consistency checks
!      do i = 1,n_cluster
!         if (clusterC(i) /= count(atomL(1:natom)==i)) then
!            stop 'error in analyse_cluster'
!         end if
!      end do
!      if (sum(clusterC(1:n_cluster))/=natom) then
!         stop 'error in analyse_cluster: wrong count'
!      end if
!   END SUBROUTINE analyse_cluster
!
!   SUBROUTINE print_cluster(natom,iu)
!      implicit none
!      integer,intent(in):: natom,iu
!      integer:: i
!      do i = 1,n_cluster
!         write(iu,*)'cluster = ',i,'count = ',clusterC(i)
!         write(iu,'(10i6)') cluster(i,1:clusterC(i))
!         if (clusterC(i) /= count(atomL(1:natom)==i)) then
!            stop 'error in print_cluster'
!         end if
!      end do
!   END SUBROUTINE print_cluster

!   SUBROUTINE delete_cluster(max_size)
!      USE connectivity_mod
!      IMPLICIT NONE
!      integer,intent(in):: max_size
!      integer:: i,j
!      integer,allocatable:: cluster_rev(:)
!
!      cluster_loop: do i = 1,n_cluster
!         if (count(cluster(i,:)/=0)<max_size) then
!            if (allocated(cluster_rev)) deallocate(cluster_rev)
!            allocate( cluster_rev(clusterC(i)) )
!            cluster_rev = 0
!            do j=1,clusterC(i)
!               cluster_rev(j) = cluster(i,clusterC(i)+1-j)
!            end do
!            do j=1,clusterC(i)
!               call delete_atom(cluster_rev(j))
!            end do
!            EXIT cluster_loop
!        end if
!      end do cluster_loop
!   END SUBROUTINE delete_cluster

!   logical function cluster_check(max_size)
!      integer,intent(in):: max_size
!      integer:: i,nc
!      nc = 0
!      do i = 1,n_cluster
!         if (clusterC(i)<max_size) then
!            nc = nc + 1
!         end if
!      end do
!      cluster_check = (nc==0)
!   end function

END MODULE HKNonLattice_mod

!!>include 'readline_mod.f90'
MODULE readline_mod
    IMPLICIT NONE
    save
    integer,parameter :: maxlnt = 132
    character(len=*),parameter :: fmtt = '(a132)'
    character(len=maxlnt) :: line
    integer :: istart,lnt
    public:: getword,getline,line,istart,lnt
!
CONTAINS

  SUBROUTINE getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
  END SUBROUTINE getline

  SUBROUTINE getword(word,ios)
      character(len=*),intent(out) :: word
      integer,intent(out):: ios
      integer :: i,iend,i1
      ios = 0
      word(:)=''
      if(istart > lnt)then
         ios = -1
         return
      end if
      i1 = 0
      do i = istart,lnt
         if(line(i:i) /= ' ')then
            i1 = i
            exit
         end if
      end do
      if(i1 == 0) then
         ios = -2
         return
      end if
      iend = 1
      do i = i1+1,lnt
         if(line(i:i) == ' ' .or. i == lnt)then
            iend = i
            exit
         end if
      end do
      if (i1 >= lnt) iend = i1
      word(1:iend-i1+1) = line(i1:iend)
      istart = i+1
      return
  END SUBROUTINE getword
END MODULE readline_mod

!!>include 'matvec3d_mod.f90'
MODULE matvec3d_mod
   USE precision_mod, only: wp
   implicit none
   public :: dotp_3d,len_3d,crossp_3d
   public :: transpose_3d,matvec_3d,matmul_3d,getinv3d
!
CONTAINS

   PURE FUNCTION dotp_3d(v1,v2)
!-----vector dot product
      real(wp):: dotp_3d
      real(wp),intent(in):: v1(3),v2(3)
      dotp_3d = v1(1)*v2(1)+ v1(2)*v2(2)+ v1(3)*v2(3)
   END FUNCTION

   PURE FUNCTION len_3d(v)
!-----absolute value of a vector
      real(wp):: len_3d
      real(wp),intent(in):: v(3)
      len_3d = sqrt( v(1)*v(1)+ v(2)*v(2)+ v(3)*v(3) )
   END FUNCTION

   PURE FUNCTION crossp_3d(v1,v2)
!-----cross product
      real(wp),dimension(3):: crossp_3d
      real(wp),dimension(3),intent(in):: v1,v2
      crossp_3d(1) = v1(2)*v2(3)-v1(3)*v2(2)
      crossp_3d(2) = v1(3)*v2(1)-v1(1)*v2(3)
      crossp_3d(3) = v1(1)*v2(2)-v1(2)*v2(1)
   END FUNCTION crossp_3d

   PURE FUNCTION transpose_3d(m)
      real(wp) :: transpose_3d(3,3)
      real(wp),intent(in) :: m(3,3)
      integer :: i,j
      forall(i=1:3,j=1:3) transpose_3d(j,i) = m(i,j)
   END FUNCTION transpose_3d

   PURE FUNCTION matvec_3d(m,v)
      real(wp) :: matvec_3d(3)
      real(wp),intent(in) :: m(3,3)
      real(wp),intent(in) :: v(3)
      matvec_3d(1) = m(1,1)*v(1) + m(1,2)*v(2) + m(1,3)*v(3)
      matvec_3d(2) = m(2,1)*v(1) + m(2,2)*v(2) + m(2,3)*v(3)
      matvec_3d(3) = m(3,1)*v(1) + m(3,2)*v(2) + m(3,3)*v(3)
   END FUNCTION matvec_3d

   PURE FUNCTION matmul_3d(m1,m2)
      real(wp),dimension(3,3) :: matmul_3d
      real(wp),dimension(3,3),intent(in) :: m1,m2
      matmul_3d(1,1) = m1(1,1)*m2(1,1) + m1(1,2)*m2(2,1) + m1(1,3)*m2(3,1)
      matmul_3d(1,2) = m1(1,1)*m2(1,2) + m1(1,2)*m2(2,2) + m1(1,3)*m2(3,2)
      matmul_3d(1,3) = m1(1,1)*m2(1,3) + m1(1,2)*m2(2,3) + m1(1,3)*m2(3,3)
      matmul_3d(2,1) = m1(2,1)*m2(1,1) + m1(2,2)*m2(2,1) + m1(2,3)*m2(3,1)
      matmul_3d(2,2) = m1(2,1)*m2(1,2) + m1(2,2)*m2(2,2) + m1(2,3)*m2(3,2)
      matmul_3d(2,3) = m1(2,1)*m2(1,3) + m1(2,2)*m2(2,3) + m1(2,3)*m2(3,3)
      matmul_3d(3,1) = m1(3,1)*m2(1,1) + m1(3,2)*m2(2,1) + m1(3,3)*m2(3,1)
      matmul_3d(3,2) = m1(3,1)*m2(1,2) + m1(3,2)*m2(2,2) + m1(3,3)*m2(3,2)
      matmul_3d(3,3) = m1(3,1)*m2(1,3) + m1(3,2)*m2(2,3) + m1(3,3)*m2(3,3)
   END FUNCTION matmul_3d

   PURE SUBROUTINE getinv3d(hmat,hmati)
      real(wp),intent(in) :: hmat(3,3)
      real(wp),intent(out) :: hmati(3,3)
      real(wp) :: det,odet
      det = hmat(1,1)*(hmat(2,2)*hmat(3,3) - hmat(2,3)*hmat(3,2)) &
          + hmat(1,2)*(hmat(2,3)*hmat(3,1) - hmat(2,1)*hmat(3,3)) &
          + hmat(1,3)*(hmat(2,1)*hmat(3,2) - hmat(2,2)*hmat(3,1))
      odet = 1.0_wp/det
      hmati(1,1) = (hmat(2,2)*hmat(3,3) - hmat(2,3)*hmat(3,2))*odet
      hmati(2,2) = (hmat(1,1)*hmat(3,3) - hmat(1,3)*hmat(3,1))*odet
      hmati(3,3) = (hmat(1,1)*hmat(2,2) - hmat(1,2)*hmat(2,1))*odet
      hmati(1,2) = (hmat(1,3)*hmat(3,2) - hmat(1,2)*hmat(3,3))*odet
      hmati(2,1) = (hmat(3,1)*hmat(2,3) - hmat(2,1)*hmat(3,3))*odet
      hmati(1,3) = (hmat(1,2)*hmat(2,3) - hmat(1,3)*hmat(2,2))*odet
      hmati(3,1) = (hmat(2,1)*hmat(3,2) - hmat(3,1)*hmat(2,2))*odet
      hmati(2,3) = (hmat(1,3)*hmat(2,1) - hmat(2,3)*hmat(1,1))*odet
      hmati(3,2) = (hmat(3,1)*hmat(1,2) - hmat(3,2)*hmat(1,1))*odet
   END SUBROUTINE getinv3d

END MODULE matvec3d_mod

!!>include 'pdb_read_mod.f90'
MODULE pdb_read_mod
   USE precision_mod
   USE global_vars_mod
   USE files_mod
   USE atom_types_mod
   USE coordinates_mod
   USE connectivity_mod
   IMPLICIT NONE
   SAVE
   integer,parameter :: maxlnt = 132
   character(len=*),parameter :: fmtt = '(a132)'
   character(len=maxlnt) :: line
   integer :: istart,lnt
   public:: getword,getline,line,istart,lnt
CONTAINS

   SUBROUTINE pdb_read(iu,natom)
      integer:: iargc,iatom,j
      integer:: ia,ig,iu,io,i,ios
      integer,intent(out):: natom
      character(len=132):: fin0,fout0,fout1
      character(len=80):: word,carg
      character(len=32):: c5,atomname
      natom = 0
      iatom = 0
      call getline(iu,ios)
      call getword(carg,ios)
      call getword(carg,ios)
      read(unit=carg,fmt=*) natom
      call getline(iu,ios)
      call getword(carg,ios)
      call getword(carg,ios)
      read(unit=carg,fmt=*) boxl(1)
      boxl(1:3) = boxl(1)
      allocate(rxyz(natom,3),proximity(natom,4),atom(natom))
      rxyz = 0.0_wp
      proximity = 0
      i=0
      lines: do
         i=i+1
         call getline(iu,ios)
         if(ios /= 0) exit
         call getword(word,ios)
         word = line(1:6)
         select case(trim(word))
         case('HETATM')
            if (i<10000) call getword(carg,ios)
            carg = line(6+1:6+5)
            read(unit=carg,fmt=*)ia
            call getword(carg,ios)
!            read(unit=carg,fmt=*)atom(name2atom(atom_name(i))
            read(unit=carg,fmt=*)atomname
            atom(i) = name2atom(trim(atomname))
            call getword(carg,ios)
            read(unit=carg,fmt=*)ig
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,1)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,2)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,3)
        case('CONECT')
            call getword(carg,ios)
            c5 = line(6+1:6+5)
            read(unit=c5,fmt=*)iatom
            do j=1,ncmax(atom(iatom))
               c5 = line(6+5*(j)+1:6+5*(j)+5)
               read(unit=c5,fmt=*) proximity(iatom,j)
            end do
            if(iatom==natom) EXIT lines
        case('END')
            EXIT
        case default
            write(*,*) 'error : ',trim(word),' is not on list'
        end select
      end do lines
   CLOSE(iu)

   END SUBROUTINE pdb_read

   SUBROUTINE getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
   END SUBROUTINE getline

   SUBROUTINE getword(word,ios)
      character(len=*),intent(out) :: word
      integer,intent(out):: ios
      integer :: i,iend,i1
      ios = 0
      word(:)=''
      if(istart > lnt)then
         ios = -1
         return
      end if
      i1 = 0
      do i = istart,lnt
         if(line(i:i) /= '')then
            i1 = i
            exit
         end if
      end do
      if(i1 == 0) then
         ios = -2
         return
      end if
      iend = 1
      do i = i1+1,lnt
         if(line(i:i) == '' .or. i == lnt)then
            iend = i
            exit
         end if
      end do
      if (i1 >= lnt) iend = i1
      word(1:iend-i1+1) = line(i1:iend)
      istart = i+1
      return
   END SUBROUTINE  getword

END MODULE pdb_read_mod

!!>include 'bond_list_mod.f90'

MODULE bond_list_mod
   use precision_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   integer:: nbondtot,nang
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars_mod, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types_mod
!      USE Keating_parameters
      integer:: ia1,ia2,j
      nbondtot = 0
      if (.not.allocated(nbond)) then
         allocate( ibond(2,natom_max*4),nbond(natom_max) )
      end if
      do ia1 = 1,natom
         nbond(ia1) = 0
         do j = 1,ncmax(atom(ia1))
            ia2 = proximity(ia1,j)
            if (ia2 > 0) nbond(ia1) = nbond(ia1) + 1
            if (ia2 > ia1) then
               nbondtot = nbondtot + 1
               ibond(1,nbondtot) = ia1
               ibond(2,nbondtot) = ia2
            end if
         end do
      end do
   END SUBROUTINE bond_list

!   subroutine set_bond_angle_lists()
!      USE global_vars_mod, only: natom,natom_max
!      USE connectivity_mod, only: proximity
!      USE atom_types_mod
!      USE Keating_parameters
!      integer:: i,ia1,ia2,ia3,ib,ic
!!
!      if (.not.allocated(ibond)) then
!         allocate( ibond(2,natom_max*4) )
!         allocate( iang(3,natom_max*6) )
!      end if
!      allocate(kang(6*natom_max),ctheta(6*natom_max))
!      nbondtot = 0
!      nang = 0
!      do ia1 = 1,natom
!      do ib = 1,ncmax(atom(ia1))
!         ia2 = proximity(ia1,ib)
!         if(ia2 == 0)cycle
!         if(ia1 < ia2)then
!            nbondtot = nbondtot+1
!            ibond(1,nbondtot) = ia1
!            ibond(2,nbondtot) = ia2
!         end if
!         do ic = 1,ncmax(atom(ia2))
!            ia3 = proximity(ia2,ic)
!            if(ia3 == 0)cycle
!            if(ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
!            if(ia1 < ia3)then
!               nang = nang+1
!               iang(1,nang) = ia1
!               iang(2,nang) = ia2
!               iang(3,nang) = ia3
!            end if
!         end do
!      end do
!      end do
!      do i = 1, nang
!         if (atom(iang(2,i)) == iSilicon) then
!            kang(i) = KOSiO
!            ctheta(i) = -1.0_wp/3.0_wp
!         else
!            kang(i) = KSiOSi
!            ctheta(i) = -1.0_wp
!         end if
!      end do
!   end subroutine set_bond_angle_lists
END MODULE bond_list_mod

!!>include 'crossbond_mod.f90'
MODULE crossbond_mod

!    READ IN THE INITIAL DATA
!    WHEN HAVE THE DATA THEN WE NEED TO IDENTIFY ALL ATOMS INVOLVED IN CROSSBONDS
!    THIS WILL BE DONE WITH THE CROSSBOND CODE IN DEVELOPEMENT. IT IS NOW NOT
!    REQUIRED TO ACTUALLY SEGREGATE THE BONDS INTO THE DIFFERENT DIRECTIONS AS THE
!    PROXIMITY VALUES OF EACH WILL NOW BE SET TO ZERO AND CALCULATED USING THE
!    NEIGHBOUR LIST ALGORITHM.
!    ALTERNATIVELY TO THE ORIGINAL APPROACH TO THE PROBLEM OF CORRECTING THE
!    CROSSBONDS TO THEIR IMAGES IN THE VARIOUS DIRECTIONS, THIS WILL IMAGE ALL ATOM
!    COORDINATES AND ALL RESPECTIVE PROXIMITY VALUES EXCEPT THOSE OF ATOMS INVOLVED
!    IN CROSSBONDS. IN THE LATTER CASE THE PROXIMITY VALUES OF THE ATOM ARE SET TO
!    ZERO.
!    FOLLOWING THIS THE SITUATION EXISTS THAT ALL THE CONNECTING ATOMS HAVE NO
!    CONNECTIVITY. THESE ARE LOOPED OVER AND AN ATTEMPT IS MADE TO CORRECT THE
!    CONNECTIVITY OF THE SYSTEM. THIS IS DONE USING THE INITIAL METHOD OF ATOMS
!    BEING WITHIN AN SI-O BOND LENGTH AWAY FROM EACH OTHER. IF IT HAPPENS THAT ONE
!    OF THESE ATOMS HAS A NUMBER OF NEIGHBOURS CLOSER THAN THE CUT-OFF DISTANCE
!    WHICH IS HIGHER THAN THE DESIRED COORDINATION THEN WE GO BACK TO THE SAME ATOM
!    IN ORIGINAL BOX AND TAKE THE "EXACT" DISTANCES BETWEEN IT AND IT'S NEIGHBOURS.
!    USING THESE DISTANCES WE SEARCH FOR NEIGHBOURS WITH CORRESPONDING EXACT
!    DISTANCES FROM TARGET ATOM.


!    if imaging a box, then it is required that the crossbond mod should be
!    initialised before the imaging takes place. then the crossbond correction
!    subroutine should be called after the imaging of the coordinates and
!    connectivity takes place.
!    NOTE:- the initial box must be a cube for this to work

   USE precision_mod
   USE coordinates_mod
   USE connectivity_mod
   IMPLICIT NONE
   integer:: ncross
   integer,allocatable:: proximity_temp(:,:)
   integer,allocatable:: cross_atom(:)
CONTAINS

   SUBROUTINE init_crossbond(n_image,boxl)
      USE bond_list_mod
      USE nlist_mod
      integer,intent(in):: n_image !total number of repeated original cells
      real(wp),intent(in):: boxl(3)

      integer:: i,j,n_crossbond
      integer,allocatable:: crossbond(:),cross_atom(:)
      real(wp):: r3(3)

      ! this will set the various values of initial box length and the proximity
      ! values required for the initial deletion of the values
      ! detail on the input box's size
      boxl2 = boxl/2.0_wp
      ! CHECK FOR CROSSBONDS.
      call bond_list
      if (.not.allocated(crossbond)) allocate(crossbond(nbondtot))
      crossbond = 0
      n_crossbond = 0
      do i=1,nbondtot ! nbondtot calc'ed in bond_list
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         r3 = SQRT(dot_product(r3,r3))
         if (  r3(1) > boxl2(1) .OR. r3(2) > boxl2(2) .OR. &
               r3(3) > boxl2(3)  ) then
            n_crossbond = n_crossbond + 1
            crossbond(n_crossbond) = i
         end if
      end do
!    store all of the atom numbers that are involved in crossbonds and the
!    connectivity of these atoms
      allocate(cross_atom((n_image)*nbondtot*2))
      j=1
      ncross = 0
      do i=1,n_crossbond
         cross_atom(j) = ibond(1,crossbond(i))
         cross_atom(j+1) = ibond(2,crossbond(i))
         ncross = ncross + 2
         j=j+2
      end do

      proximity_temp(1:natom,:) = proximity(1:natom,:)

      do i =1, ncross
         proximity(cross_atom(i),:) = 0
      end do
   END SUBROUTINE init_crossbond


   SUBROUTINE crossbond_correct(nix,niy,niz,boxl_0,natom_orig)
      !!! NOTE:- only works for atoms with 4 or less neighbours 
      USE precision_mod, only: wp
      USE global_vars_mod
      USE bond_list_mod
      USE atom_types_mod
      USE connectivity_mod
      USE files_mod
      USE nlist_mod
      ! nix, niy, & niz are the number of images in each respective direction so
      ! need to have this passed in...probably the easiest thing to do
      integer,intent(in):: natom_orig
      real(wp),intent(in):: nix,niy,niz,boxl_0(3)

      integer:: ic,nc,iatom
      integer:: i,ii,k,kk,m,n !counters and dummy args
      integer:: ifirst
      real(wp):: r3(3),bond_lens(4),dr(3)
      logical:: consistent

      !  need to specify the changes in the box length directions for
      !  initialisation of the neighbour list.
      boxl(1:3) = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
      boxl2(1:3) = boxl(1:3)/2.0_wp

      if (allocated(ll_old)) deallocate(ll_old)
      if (allocated(ll)) deallocate(ll)
      if (allocated(lr)) deallocate(lr)
      if (allocated(HOC)) deallocate(HOC)
      if (allocated(HOC_OLD)) deallocate(HOC_OLD)

      print*, 'nix = ',nix
      print*, 'niy = ',niy
      print*, 'niz = ',niz

      call INIT_NLIST(boxl,5.0_wp)
      call NEW_NLIST

      ! the following is to get the connectivity of the edge atoms using the
      ! adjusted nlist
      main_loop: do i = 1, natom
         if ( .NOT. ALL( proximity(i,:)==0 ) ) cycle
         r3(1:3) = rxyz(i,1:3)
         ic = CELL(r3)
         CALL NEIGCELL(ic,2,neigh,ncell)
         n=0
         if (i>natom_orig) then
            iatom = mod(i,natom_orig)
         elseif (i<=natom_orig) then
            iatom = i
         end if
         ! need to have an array containing the connectivity
         ! of the target atom i. iatom is the original of
         ! the image under investigation.
         boxl = (/boxl_0(1),boxl_0(2),boxl_0(3)/)
         boxl2 = boxl/2.0_wp
         bond_lens = 0.0_wp
         do k = 1, ncmax(atom(iatom))
            dr = rxyz(iatom,:) - rxyz(proximity_temp(iatom,k),:)
            CALL pbc(dr)
            bond_lens(k) = dot_product(dr,dr)
         end do
         !reset boxl to proper value
         boxl = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
         boxl2 = boxl/2.0_wp
            inner_cell_loop: do kk=1,neigh
               nc=ncell(kk)
               if (nc==0) cycle inner_cell_loop
               m=HOC(nc)
               inner_cell_atom_loop: do while (m/=0)
                  dr = r3-rxyz(m,1:3)
                  if (m==i) GOTO 20
                  CALL pbc(dr)
                  do ii = 1, ncmax(atom(i))
                     if (ABS( dot_product(dr,dr)-bond_lens(ii) )<=0.00001_wp) EXIT
                     if (ii==ncmax(atom(i))) GOTO 20
                  end do
                  n=n+1
                  proximity(i,n) = m
                  if (n>4) stop 'need better condition'
                  if (n==ncmax(atom(i))) cycle main_loop
20                CONTINUE
                  m=LL(m)
               end do inner_cell_atom_loop
            end do inner_cell_loop
      end do main_loop

      boxl = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
      boxl2 = boxl/2.0_wp
      call bond_list

! turn on the following checks if you don't trust the correction
!      call check_proximity(consistent,ifirst,k)
!      if (.NOT. consistent) then
!         print*, ''
!         print*, 'not consistent'
!         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
!         print*, 'proximity(',k,',:)=', proximity(k,:)
!      end if
!
!      call check_proximity2(consistent,ifirst,k)
!      if (.NOT. consistent) then
!         print*, ''
!         print*, 'not consistent_2'
!         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
!         print*, 'atom_name(atom(ifirst)) = ', atom_name(atom(ifirst))
!         print*, 'atom_name(atom(k)) = ', atom_name(atom(k))
!      end if
   END SUBROUTINE crossbond_correct
END MODULE crossbond_mod


!!>include 'command_line_mod.f90'
MODULE COMMAND_LINE_MOD
   implicit none
CONTAINS

   function command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   end function

   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
      CHARACTER(LEN=1000) :: TMPVAL
      INTEGER,EXTERNAL :: IARGC
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      ELSE IF (NUMBER > IARGC()) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 2
          RETURN
      END IF
      IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
      IF (PRESENT(LENGTH)) THEN
          IF (PRESENT(VALUE)) THEN
              LENGTH = LEN_TRIM(VALUE)
          ELSE
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          END IF
      END IF
      IF (PRESENT(STATUS)) STATUS = 0
      RETURN
   END SUBROUTINE GET_COMMAND_ARGUMENT

END MODULE COMMAND_LINE_MOD




PROGRAM etch_amorphous
      USE precision_mod
      USE atom_types_mod
      USE files_mod
      USE coordinates_mod
      USE connectivity_mod
      USE global_vars_mod
      USE nlist_mod
      USE HKNonlattice_mod
      USE readline_mod
      USE sort_mod
      USE matvec3d_mod
      USE pdb_read_mod
      USE crossbond_mod
      USE COMMAND_LINE_MOD
   
      IMPLICIT NONE
      integer,parameter:: iBorate=0
      integer,parameter:: iSilicate=1
      real(wp),parameter:: bondl_SiO = 1.60_wp
      character(2),parameter:: atom_name2(0:1)=(/ 'B ','Si' /)
   
      integer:: narg,leng,stat
      integer:: i,j,k,ii,jj,iu,ic,nc
      integer:: floater
      integer:: ia,ig,ios,iatom
      integer:: ps_natom,ndel,tmp
      integer:: n_OH,coord_type,nat_temp
      integer:: ifirst,i2,ucoord1,ucoord2,ucoord3,MaxCluster
      integer:: nSi_1_OH,nSi_1_OH2,nSi_1_OH3
      integer:: nSi_2_OH,nSi_2_OH2,nSi_2_OH3
      integer:: nSi_3_OH,nSi_3_OH2,nSi_3_OH3
      integer:: nSi_4_OH,nSi_4_OH2,nSi_4_OH3
      integer:: nimage,natom_orig
      integer,allocatable:: ps_atom(:),del_list(:),ClusterCount(:)
   
      real(wp):: x,y,z,dr(3),nix,niy,niz,boxl_0(3)
      real(wp):: cutoff,cutoff2
      real(wp):: t(0:10)
      real(wp):: ra(3),rb(3),rc(3),r_ba(3),r_ca(3),mag_crossp_3d,temp(3)
      real(wp),allocatable:: ps_rxyz(:,:)
   
      logical:: consistent
      logical,allocatable:: dl_logical(:)
   
      character(len=132):: finput1,finput2,finput3,finput4,carg
      character(len=32):: ctmp,word,atomname,c5
   
      nseed=0
   
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      if(narg /= 2) then
         write(*,*)'usage :'
         write(*,'(/)')
         write(*,*)'  executable_name phase_sep.xyz amorphous.xyz amorphous.conn'
         write(*,'(/)')
      end if
   
      call get_command_argument (1, finput1, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
      call get_command_argument (2, finput2, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
   
      call get_command_argument (3, finput3, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
   
   
      cutoff = 6.0_wp
      cutoff2 = cutoff**2
   
      ! read in the phase separation data
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput1),STATUS='old')
         read(iu,*) ps_natom
         allocate(ps_rxyz(ps_natom,3),ps_atom(ps_natom))
         ps_rxyz = 0.0_wp
   
         read(iu,*) ctmp
         do i = 1,ps_natom
            read(iu,*) ctmp,x,y,z
            ps_atom(i) = name2atom2(ctmp)
            ps_rxyz(i,1:3) = (/ x,y,z /)
         end do
      CLOSE(iu)
   
   
      ! read in coordinates of file amorphous silica system
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput2),STATUS='old')
         read(iu,*) natom
         allocate(rxyz(2*natom,3),proximity(2*natom,4),&
                  proximity_temp(2*natom,4),atom(2*natom))
         natom_orig=natom
         proximity=0

         read(iu,*) boxl(1),boxl(2),boxl(3)
         do i = 1,natom
            read(iu,*) ctmp,x,y,z
            atom(i)=name2atom(ctmp)
            rxyz(i,1:3) = (/ x,y,z /)
         end do
      CLOSE(iu)
   
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput3),STATUS='old')
         call read_conn(iu)
      CLOSE(iu)
   
      natom_max = 2*natom
   
      call prox_org
      
      
      call check_proximity2(consistent,ifirst,i2)
      if (.NOT. consistent) then
      print'(2/)'
      print*,'pre-crossbond'
      print*, 'not consistent-wrong atoms bonded'
      print'(/)'
      print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
      print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
      print'(/)'
      print*, 'atom(',i2,') = ',atom_name(atom(i2))
      print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
      print'(3/)'
      stop
      end if
      
      call check_proximity(consistent,ifirst,i2)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent'
      print'(/)'
      print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
      print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
      print'(/)'
      print*, 'atom(',i2,') = ',atom_name(atom(i2))
      print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
      print'(3/)'
      stop
      end if
   
      boxl_0(1:3)=boxl(1:3)
      boxl2(1:3)=boxl(1:3)/2.0_wp
      natom_max = 2*natom
      allocate(del_list(2*natom))
      del_list=0
      call init_nlist(boxl,5.0_wp)
      call new_nlist
      call cpu_time(t(0))
      ndel = 0
      del_list = 0
      allocate(dl_logical(natom_max))
      dl_logical(1:natom_max)=.false.
   
      ! loop which marks atoms that overlap the 'boron' phase as .TRUE.
      ! for delete
      atom_loop: do i = 1,ps_natom
         !if (mod(i,100)==0) print*,'i = ', i
         if (ps_atom(i) == iBorate) then
            ic = CELL(ps_rxyz(i,1:3))
            call NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj=1,neigh
               nc = ncell(jj)
               if (nc == 0) cycle cell_loop
               j = HOC(nc)
               cell_atom_loop: do while (j /= 0)
                  dr(1:3) = ps_rxyz(i,1:3) - rxyz(j,1:3)
                  call pbc(dr)
                  if ( dot_product(dr,dr) < cutoff2 )  then
                     dl_logical(j)=.TRUE.
                  end if
1020              j = LL(j)
               end do cell_atom_loop
            end do cell_loop
         end if
      end do atom_loop
   
   
   
      ! do a simple imaging of the cube in the negative z direction
      ! first move initial centred on origin up a length boxl2(3)
      nimage=2
      rxyz(1:natom,3) = rxyz(1:natom,3) + boxl2(3)
      CALL init_crossbond(nimage,boxl)
      rxyz(natom+1:2*natom,1) = rxyz(1:natom,1)
      rxyz(natom+1:2*natom,2) = rxyz(1:natom,2)
      rxyz(natom+1:2*natom,3) = rxyz(1:natom,3) - boxl(3)
      atom(natom+1:2*natom) = atom(1:natom)
      where (proximity(1:natom,:)/=0) proximity(natom+1:2*natom,:) = &
                      proximity(1:natom,:) + natom
   
      ncross=2*ncross
      natom=2*natom
      nix=1.0_wp
      niy=1.0_wp
      niz=2.0_wp
      CALL crossbond_correct(nix,niy,niz,boxl_0,natom_orig)
   
      call prox_org
   
      ! a check to make sure the imaging and thus connectivity corrections
      ! were done correctly   
      call check_proximity(consistent,ifirst,i2)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent'
      print'(/)'
      print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
      print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
      print'(/)'
      print*, 'atom(',i2,') = ',atom_name(atom(i2))
      print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
      print'(3/)'
      stop
      end if
      
      call check_proximity2(consistent,ifirst,i2)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent-wrong atoms bonded - post crossbond'
      print'(/)'
      print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
      print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
      print'(/)'
      print*, 'atom(',i2,') = ',atom_name(atom(i2))
      print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
      print'(3/)'
      stop
      end if

      ! delete atoms
      print*,' ndel = ', count(dl_logical(1:natom)==.TRUE.)
      do i = natom,1,-1
         if (dl_logical(i)) call delete_atom(i)
      end do

call cpu_time(t(1))
print*, 'delete time (s): ', t(1)-t(0)

   
      ! remove floating atoms
      del_list=0
      floater = 0
      do k = 1, natom
         if (count(proximity(k,:)/=0) == 0) then
            floater = floater + 1
            del_list(floater)=k
         end if
      end do
      call shell(floater,del_list)
   
     ! sort and delete atoms
      call shell(floater,del_list)
      do i = floater,1,-1
         call delete_atom(del_list(i))
      end do
   
print*, 'floater = ', floater
print*, 'natom pre-float del',natom

      call re_init_nlist(boxl,5.0_wp)
      call new_nlist

      !! count number of silicons which are not 4 coordinated
      !! and which have no OH group
      ucoord1 = 0
      ucoord2 = 0
      ucoord3 = 0
      coord_loop: do i = 1,natom
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) > 0) then
            do j=1,4
               if (proximity(i,j)/=0) then
                  if (atom(proximity(i,j))==iOxygenH) CYCLE coord_loop
               end if
            end do
            select case(count(proximity(i,:)/=0))
            case(1)
               ucoord1 = ucoord1 + 1
            case(2)
               ucoord2 = ucoord2 + 1
            case(3)
               ucoord3 = ucoord3 + 1
            end select
         end if
         end if
      end do coord_loop
      print*, 'ucoord (Si_0_1) = ', ucoord1
      print*, 'ucoord (Si_0_2) = ', ucoord2
      print*, 'ucoord (Si_0_3) = ', ucoord3
      print'(2/)'
   
      call prox_org

do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                                      !
      !              remove any 1 coord Si                   !
      !              remove any 2 coord Si                   !
      !           add OH to any 3 coord Si                   !
      !                                                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ndel = 0
      del_list = 0
      nat_temp = natom
      at_loop: do i = 1,nat_temp
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) > 0) then
   
            select case(count(proximity(i,:)/=0))
            case(1)
               do j=1,ncmax(atom(i))
                  if (proximity(i,j)==0) cycle
                  if (atom(proximity(i,j))==iOxygen) then
                     atom(proximity(i,j)) = iOxygenH
                  end if
               end do
               ndel = ndel + 1
               del_list(ndel) = i
   
            case(2)
               do j=1,ncmax(atom(i))
                  if (proximity(i,j)==0) cycle
                  if (atom(proximity(i,j))==iOxygen) then
                     atom(proximity(i,j)) = iOxygenH
                  end if
               end do
               ndel = ndel + 1
               del_list(ndel) = i
   
            case(3)
               ! add on an OH to the empty spot
               natom = natom + 1
               ra(1:3) = rxyz(proximity(i,1),:)
               rb(1:3) = rxyz(proximity(i,2),:)
               rc(1:3) = rxyz(proximity(i,3),:)
               r_ba = ra - rb
               call pbc(r_ba)
               r_ca = ra - rc
               call pbc(r_ca)
               mag_crossp_3d = sqrt(dot_product(crossp_3d(r_ba,r_ca),crossp_3d(r_ba,r_ca)))
               temp = crossp_3d(r_ba,r_ca)/mag_crossp_3d
               rxyz(natom,:) = rxyz(i,:) + temp*bondl_SiO
               call pbc(rxyz(natom,:))
               atom(natom) = iOxygenH
               proximity(natom,1) = i
               proximity(i,4) = natom
            end select
         end if
         end if
      end do at_loop
   
      call shell(ndel,del_list)
   
print*, 'ndel = ', ndel
print*, 'natom pre-treat - del',natom

      do k=ndel,1,-1
         call delete_atom(del_list(k))
      end do
print*, 'natom post-treat - del',natom
print'(2/)'

      call prox_org
   
      ! remove floating atoms
      del_list=0
      floater = 0
      do k = 1, natom
         if (count(proximity(k,:)/=0) == 0) then
            floater = floater + 1
            del_list(floater)=k
         end if
      end do
      call shell(floater,del_list)
   
     ! sort and delete atoms
      call shell(floater,del_list)
      do i = floater,1,-1
         call delete_atom(del_list(i))
      end do

print*, 'floater = ', floater

      CALL Init_HKNonLattice(natom)
      CALL HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
      del_list=0
      ndel=0
      do i=1,natom
         if (atomL(i)/=MaxCluster) then
            ndel=ndel+1
            del_list(ndel)=i
         end if
      end do

      call shell(ndel,del_list)
      do i=ndel,1,-1
         call delete_atom(del_list(i))
      end do
      deallocate(ClusterCount)


      !! count number of silicons which are not 4 coordinated
      !! and which have no OH group
      ucoord1 = 0
      ucoord2 = 0
      ucoord3 = 0
      do i = 1,natom
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) > 0) then
!            if (ANY(atom(proximity(i,:))==iOxygenH)) CYCLE
            do j=1,4
               if (proximity(i,j)/=0) then
                  if (atom(proximity(i,j))==iOxygenH) CYCLE
               end if
            end do
   
            select case(count(proximity(i,:)/=0))
            case(1)
               ucoord1 = ucoord1 + 1
            case(2)
               ucoord2 = ucoord2 + 1
            case(3)
               ucoord3 = ucoord3 + 1
            end select
         end if
         end if
      end do
print*, 'ucoord (Si_0_1) = ', ucoord1
print*, 'ucoord (Si_0_2) = ', ucoord2
print*, 'ucoord (Si_0_3) = ', ucoord3
print'(2/)'

      if (ucoord1==0 .AND. ucoord2==0 .AND. ucoord3==0) EXIT
end do

      del_list=0
      floater = 0
      do k = 1, natom
         if (count(proximity(k,:)/=0) == 0) then
            floater = floater + 1
            del_list(floater)=k
         end if
      end do
      call shell(floater,del_list)

print*, 'floater = ', floater
print*, 'natom pre-float del',natom

      do k=floater,1,-1
         call delete_atom(del_list(k))
      end do
print*, 'natom post-float del',natom

     call prox_org
 
     ! go over and delete any clusters of floating atoms

      CALL Init_HKNonLattice(natom)
      CALL HKNonLattice(natom,proximity,n_cluster,atomL)
      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
      del_list=0
      ndel=0
      do i=natom,1-1
         if (atomL(i)/=MaxCluster) then
            ndel=ndel+1
            del_list(ndel)=i
         end if
      end do

      do i=1,ndel
         call delete_atom(del_list(i))
      end do

do
      !! delete any FULLY bonded si atoms which have an 3 OH's
      ndel = 0
      del_list = 0
      do i = 1,natom
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) == 0) then
            select case(count(atom(proximity(i,:))==iOxygenH))
            case(1)
               ! do nothing
            case(2)
               ! do nothing
            case(3)
            ndel = ndel + 1
            del_list(ndel) = i
   
            do j = 1,ncmax(atom(i))
               if (proximity(i,j)==0) stop 'screwed up'
               if (atom(proximity(i,j)) == iOxygenH) then
                  ndel = ndel + 1
                  del_list(ndel) = proximity(i,j)
               else if (atom(proximity(i,j)) == iOxygen) then
                  atom(proximity(i,j)) = iOxygenH
               end if
            end do
            case default
               cycle
            end select
         end if
         end if
      end do
   
      call shell(ndel,del_list)
   
      do k=ndel,1,-1
        call delete_atom(del_list(k))
      end do
   
      !! delete any TRIPLE bonded si atoms which have an 3 OH's
      ndel = 0
      del_list = 0
      do i = 1,natom
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) == 1) then
            ndel = ndel + 1
            del_list(ndel) = i
            select case(count(atom(proximity(i,:))==iOxygenH))
            case(1)
               ! do nothing
            case(2)
               ! no nothing
            case(3)
               do j = 1,ncmax(atom(i))
                  if (proximity(i,j)==0) CYCLE
                  if (atom(proximity(i,j)) == iOxygenH) then
                     ndel = ndel + 1
                     del_list(ndel) = proximity(i,j)
                  else if (atom(proximity(i,j)) == iOxygen) then
                     atom(proximity(i,j)) = iOxygenH
                  end if
               end do
            end select
   
         end if
         end if
      end do
   
      call shell(ndel,del_list)
      do k=ndel,1,-1
         call delete_atom(del_list(k))
      end do
print*, 'natom post-triple - del',natom
      
      !! convert all singly bonded oxygens to
      do i = 1,natom
         if (atom(i)/=iOxygen) CYCLE
         if (count(proximity(i,:)==0) == 3) then
            atom(i) = iOxygenH
         end if
      end do

print*, 'post convert'

      !! delete any double bonded si atoms which have an OH
      ndel = 0
      del_list = 0
      do i = 1,natom
         if (atom(i) == iSilicon) then
         if (count(proximity(i,:)==0) == 2) then
            do k=1,4
               if (proximity(i,k)/=0) then
                  if (atom(proximity(i,k))==iOxygenH) then
                     ndel = ndel + 1
                     del_list(ndel) = i
                     do j = 1,ncmax(atom(i))
                        if (proximity(i,j)==0) cycle
                        if (atom(proximity(i,j)) == iOxygenH) then
                           ndel = ndel + 1
                           del_list(ndel) = proximity(i,j)
                        else if (atom(proximity(i,j)) == iOxygen) then
                           atom(proximity(i,j)) = iOxygenH
                        end if
                     end do
                  end if
               end if
            end do
         end if
         end if
      end do

      call shell(ndel,del_list)

print*, 'ndel = ', ndel
print*, 'natom double-del',natom
   
      do k=ndel,1,-1
         call delete_atom(del_list(k))
      end do
      !! delete any singly bonded si atoms
      ndel = 0
      del_list = 0
      do i = 1,natom
         if (atom(i)==iSilicon) then
         if (count(proximity(i,:)==0) == 3) then
            ndel = ndel + 1
            del_list(ndel) = i
            do j = 1,ncmax(atom(i))
               if (proximity(i,j)==0) cycle
               if (atom(proximity(i,j)) == iOxygen) then
                  atom(proximity(i,j)) = iOxygenH
               end if
            end do
         end if
         end if
      end do

      !! analyse the bonds in the system
      nSi_1_OH = 0
      nSi_1_OH2 = 0
      nSi_1_OH3 = 0
      nSi_2_OH = 0
      nSi_2_OH2 = 0
      nSi_2_OH3 = 0
      nSi_3_OH = 0
      nSi_3_OH2 = 0
      nSi_3_OH3 = 0
      nSi_4_OH = 0
      nSi_4_OH2 = 0
      nSi_4_OH3 = 0   
      do i=1,natom
         if (atom(i) == iSilicon) then
            n_OH = 0
            do j = 1,ncmax(atom(i))
               if (atom(proximity(i,j)) == iOxygenH) then
                  n_OH = n_OH + 1
               end if
            end do
            coord_type = count(proximity(i,:)/=0)
            select case(coord_type)
            case(1)
               select case(n_OH)
               case(1)
                  nSi_1_OH = nSi_1_OH + 1
               case(2)
                  nSi_1_OH2 = nSi_1_OH2 + 1
               case(3)
                  nSi_1_OH3 = nSi_1_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(2)
               select case(n_OH)
               case(1)
                  nSi_2_OH = nSi_2_OH + 1
               case(2)
                  nSi_2_OH2 = nSi_2_OH2 + 1
               case(3)
                  nSi_2_OH3 = nSi_2_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(3)
               select case(n_OH)
               case(1)
                  nSi_3_OH = nSi_3_OH + 1
               case(2)
                  nSi_3_OH2 = nSi_3_OH2 + 1
               case(3)
                  nSi_3_OH3 = nSi_3_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(4)
               select case(n_OH)
               case(1)
                  nSi_4_OH = nSi_4_OH + 1
               case(2)
                  nSi_4_OH2 = nSi_4_OH2 + 1
               case(3)
                  nSi_4_OH3 = nSi_4_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            end select
         else if (atom(i) == iOxygen) then
            ! do nothing
         else if (atom(i) == iOxygenH) then
            ! do nothing
         else
            print*, 'error in atom types'
            print*, 'atom(',i,') = ', atom_name(atom(i))
            print'(2/)'
            print*, 'program terminated prematurely'
            print'(/)'
            stop
         end if
      end do
   
      !! bonding summary
      print'(2/)'
!      print*, 'nSi_0H   = ', nSi_OH
!      print*, 'nSi_0H2  = ', nSi_OH2
!      print*, 'nSi_0H3  = ', nSi_OH3
      print*, 'Singly bonded silicons'
      print*, 'nSi_1_OH  = ',  nSi_1_OH
      print*, 'nSi_1_OH2 = ',  nSi_1_OH2
      print*, 'nSi_1_OH3 = ',  nSi_1_OH3
      
      print'(2/)'
      print*, 'Double bonded silicons'
      print*, 'nSi_2_OH  = ',  nSi_2_OH
      print*, 'nSi_2_OH2 = ',  nSi_2_OH2
      print*, 'nSi_2_OH3 = ',  nSi_2_OH3
      
      print'(2/)'
      print*, 'Triple bonded silicons'
      print*, 'nSi_3_OH  = ',  nSi_3_OH
      print*, 'nSi_3_OH2 = ',  nSi_3_OH2
      print*, 'nSi_3_OH3 = ',  nSi_3_OH3
      
      print'(2/)'
      print*, 'Fully bonded silicons'
      print*, 'nSi_4_OH  = ',  nSi_4_OH
      print*, 'nSi_4_OH2 = ',  nSi_4_OH2
      print*, 'nSi_4_OH3 = ',  nSi_4_OH3
      print'(2/)'

      ! if all ducks are in a row break out
      if (nSi_1_OH==0 .AND.&
          nSi_1_OH2==0 .AND.&
          nSi_1_OH3==0 .AND.&
          nSi_2_OH==0 .AND.&
          nSi_2_OH2==0 .AND.&
          nSi_2_OH3==0 .AND.&
          nSi_3_OH==0 .AND.&
          nSi_3_OH2==0 .AND.&
          nSi_3_OH3==0 .AND.&
          nSi_4_OH3==0) exit
end do
   
   
      !! convert all singly bonded oxygens to
      do i = 1,natom
         if (atom(i)/=iOxygen) CYCLE
         if (count(proximity(i,:)==0) == 3) then
            atom(i) = iOxygenH
         end if
      end do
   
      del_list=0
      floater = 0
      do k = 1, natom
         if (count(proximity(k,:)/=0) == 0) then
            floater = floater + 1
            del_list(floater)=k
         end if
      end do
      call shell(floater,del_list)
   
   print*, 'floater = ', floater
   print*, 'natom pre-float del',natom
   
      do k=floater,1,-1
         call delete_atom(del_list(k))
      end do
   
    
      !! analyse the bonds in the system
      nSi_1_OH = 0
      nSi_1_OH2 = 0
      nSi_1_OH3 = 0
      nSi_2_OH = 0
      nSi_2_OH2 = 0
      nSi_2_OH3 = 0
      nSi_3_OH = 0
      nSi_3_OH2 = 0
      nSi_3_OH3 = 0
      nSi_4_OH = 0
      nSi_4_OH2 = 0
      nSi_4_OH3 = 0
      do i=1,natom
         if (atom(i) == iSilicon) then
            n_OH = 0
            do j = 1,ncmax(atom(i))
               if (atom(proximity(i,j)) == iOxygenH) then
                  n_OH = n_OH + 1
               end if
            end do
      
      
            coord_type = count(proximity(i,:)/=0)
            select case(coord_type)
            case(1)
               select case(n_OH)
               case(1)
                  nSi_1_OH = nSi_1_OH + 1
               case(2)
                  nSi_1_OH2 = nSi_1_OH2 + 1
               case(3)
                  nSi_1_OH3 = nSi_1_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(2)
               select case(n_OH)
               case(1)
                  nSi_2_OH = nSi_2_OH + 1
               case(2)
                  nSi_2_OH2 = nSi_2_OH2 + 1
               case(3)
                  nSi_2_OH3 = nSi_2_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(3)
               select case(n_OH)
               case(1)
                  nSi_3_OH = nSi_3_OH + 1
               case(2)
                  nSi_3_OH2 = nSi_3_OH2 + 1
               case(3)
                  nSi_3_OH3 = nSi_3_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(4)
               select case(n_OH)
               case(1)
                  nSi_4_OH = nSi_4_OH + 1
               case(2)
                  nSi_4_OH2 = nSi_4_OH2 + 1
               case(3)
                  nSi_4_OH3 = nSi_4_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            end select
      
      
         else if (atom(i) == iOxygen) then
            ! do nothing
      
         else if (atom(i) == iOxygenH) then
            ! do nothing
      
         else
            print*, 'error in atom types'
            print*, 'atom(',i,') = ', atom_name(atom(i))
            print'(2/)'
            print*, 'program terminated prematurely'
            print'(/)'
            stop
         end if
      end do
         
      !! bonding summary
      print'(2/)'
      ! print*, 'nSi_0H   = ', nSi_OH
      ! print*, 'nSi_0H2  = ', nSi_OH2
      ! print*, 'nSi_0H3  = ', nSi_OH3
      print*, 'Singly bonded silicons'
      print*, 'nSi_1_OH  = ',  nSi_1_OH
      print*, 'nSi_1_OH2 = ',  nSi_1_OH2
      print*, 'nSi_1_OH3 = ',  nSi_1_OH3
      
      print'(2/)'
      print*, 'Double bonded silicons'
      print*, 'nSi_2_OH  = ',  nSi_2_OH
      print*, 'nSi_2_OH2 = ',  nSi_2_OH2
      print*, 'nSi_2_OH3 = ',  nSi_2_OH3
      
      print'(2/)'
      print*, 'Triple bonded silicons'
      print*, 'nSi_3_OH  = ',  nSi_3_OH
      print*, 'nSi_3_OH2 = ',  nSi_3_OH2
      print*, 'nSi_3_OH3 = ',  nSi_3_OH3
      
      print'(2/)'
      print*, 'Fully bonded silicons'
      print*, 'nSi_4_OH  = ',  nSi_4_OH
      print*, 'nSi_4_OH2 = ',  nSi_4_OH2
      print*, 'nSi_4_OH3 = ',  nSi_4_OH3
      print'(2/)'
      
      call prox_org
      
!      do i=1,natom
!      if (count(proximity(i,:)/=0)<ncmax(atom(i))) then
!      if (atom(i)==iOxygenH) then
!      if (count(proximity(i,:)/=0)==1) cycle
!      end if
!            print*, 'atom(',i,') = ', atom_name(atom(i))
!            print'(a5,1x,i6,a6,1x,4i6)', 'prox(',i,',:) = ',(proximity(i,j),j = 1,ncmax(atom(i)))
!            print'(2/)'
!         end if
!      end do
!
      call check_proximity(consistent,ifirst,i2)
      if (.NOT. consistent) then
         print'(2/)'
         print*, 'not consistent 1 - end'
         print'(/)'
         print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
         print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
         print'(/)'
         print*, 'atom(',i2,') = ',atom_name(atom(i2))
         print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
         print'(3/)'
         stop
      end if
      
      call prox_org
      
      call check_proximity2(consistent,ifirst,i2)
      if (.NOT. consistent) then
         print'(2/)'
         print*, 'not consistent 2 - end'
         print'(/)'
         print*, 'atom(',ifirst,') = ',atom_name(atom(ifirst))
         print'(a10,i7,a6,4i7)', 'proximity(',ifirst,',:) = ',proximity(ifirst,:)
         print'(/)'
         print*, 'atom(',i2,') = ',atom_name(atom(i2))
         print'(a10,i7,a6,4i7)', 'proximity(',i2,',:) = ',proximity(i2,:)
         print'(3/)'
         stop
      end if


   call get_free_file_unit(iu)
   OPEN(UNIT=iu,FILE='large_etch.xyz',STATUS='UNKNOWN')
      write(iu,*) natom
      write(iu,*) 'Etched_Structure'
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)), rxyz(i,1:3)
      end do
   CLOSE(iu)

   call get_free_file_unit(iu)
   OPEN(UNIT=iu,FILE='large_etch.conn',STATUS='UNKNOWN')
      write(iu,*) natom
      write(iu,*) 'connectivity data for amorphous.xyz'
      do i = 1,natom
         write(iu,'(5i6)') i,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
   CLOSE(iu)

110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3)



CONTAINS

   SUBROUTINE prox_org
      ! arranges the non-zero elements of each proximity
      ! array from the leftmost side
      USE connectivity_mod, only: proximity
      USE global_vars_mod, only: natom
      IMPLICIT NONE
      integer:: i,ptmp(4),nat
      do i=1,natom
         ptmp(1:4) = proximity(i,1:4)
         proximity(i,:) = 0
         nat = 0
         do j=1,4
            if (ptmp(j)/=0) then
               nat = nat + 1
               proximity(i,nat) = ptmp(j)
            end if
         end do
      end do
   END SUBROUTINE prox_org



   pure function name2atom2(c)
      integer:: name2atom2
      character(*),intent(in):: c
      integer:: i
      do i = 0,1
         if(trim(c) == trim(atom_name2(i)))then
            name2atom2 = i
            exit
         end if
      end do
      if (i > ntyp) name2atom2 = -1
   end function

   SUBROUTINE read_conn(iu)
      USE connectivity_mod
      USE atom_types_mod
      integer,intent(in):: iu
      integer:: nat,i,j,iatom
      character(len=132):: line
      character(len=80):: ctmp
      character(len=32):: c5

      ! read in a formatted connectivity file
      read(iu,*) nat
      read(iu,*) ctmp
      do i=1,nat
         line=''
         read(unit=iu,fmt='(a132)',iostat=ios) line
         c5=line(1:6)
         read(unit=c5,fmt=*)iatom
         do j=1,ncmax(atom(iatom))
            c5=line(6*(j)+1:6*(j)+6)
            read(unit=c5,fmt=*) proximity(iatom,j)
         end do
      end do

   END SUBROUTINE read_conn


END PROGRAM etch_amorphous

