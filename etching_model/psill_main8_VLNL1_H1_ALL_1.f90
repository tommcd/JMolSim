
!!>include 'precision_mod.f90'

MODULE precision_mod
   implicit none
!  integer, parameter :: sp = kind(1.0)
!  integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: sp = selected_real_kind(6,30)
   integer, parameter :: dp = selected_real_kind(15,300)
   integer, parameter :: qp_preferred = selected_real_kind(30,1000)
   integer, parameter :: qp = (1 + sign(1,qp_preferred))/2*qp_preferred+ &
                              (1 - sign(1,qp_preferred))/2*dp
!
   integer,parameter,public :: wp = dp
   integer,parameter,public :: i4b = selected_int_kind(9)
END MODULE precision_mod

!!>include 'global_vars.f90'

MODULE global_vars
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, etot
   integer:: natom,nseed,natom_max
   integer:: nattached, nrelax
END MODULE

!!>include 'sort_mod.f90'

MODULE sort_mod
      implicit none
! Quicksort modified from Numerical Recipes
CONTAINS

   SUBROUTINE qsort(n,arr0,arr)
      integer,intent(in):: n,arr0(:)
      integer,intent(out):: arr(:)
      integer,parameter:: M = 7,NSTACK = 50
      integer:: i,ir,j,jstack,k,L,istack(NSTACK)
      integer:: a,temp
      jstack = 0
      L = 1
      ir = n
      arr(1:n) = arr0(1:n)
1     if (ir - L < M) then
        do j = L + 1,ir
          a = arr(j)
          do i = j - 1,1, -1
            if (arr(i) <= a) GOTO 2
            arr(i + 1) = arr(i)
          end do
          i = 0
2         arr(i + 1) = a
        end do
        if (jstack == 0) RETURN
        ir = istack(jstack)
        L = istack(jstack - 1)
        jstack = jstack - 2
      else
        k = (L + ir)/2
        temp = arr(k)
        arr(k) = arr(L + 1)
        arr(L + 1) = temp
        if (arr(L + 1) > arr(ir)) then
          temp = arr(L + 1)
          arr(L + 1) = arr(ir)
          arr(ir) = temp
        end if
        if (arr(L) > arr(ir)) then
          temp = arr(L)
          arr(L) = arr(ir)
          arr(ir) = temp
        end if
        if (arr(L + 1) > arr(L)) then
          temp = arr(L + 1)
          arr(L + 1) = arr(L)
          arr(L) = temp
        end if
        i = L + 1
        j = ir
        a = arr(L)
3       continue
          i = i + 1
        if (arr(i) < a) GOTO 3
4       continue
          j = j - 1
        if (arr(j) > a) GOTO 4
        if (j < i) GOTO 5
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        GOTO 3
5       arr(L) = arr(j)
        arr(j) = a
        jstack = jstack + 2
        if (jstack > NSTACK) stop 'NSTACK too small in qsort'
        if (ir - i + 1 >= j - L) then
          istack(jstack) = ir
          istack(jstack - 1) = i
          ir = j - 1
        else
          istack(jstack) = j - 1
          istack(jstack - 1) = L
          L = i
        end if
      end if
      GOTO 1
   END SUBROUTINE
!  (C) Copr. 1986-92 Numerical Recipes Software #1-0zV'n26) B3.


   SUBROUTINE shell(n,v)
      integer,intent(in):: n
      integer,intent(inout):: v(:)
! Sorts vector v(1:n) into ascending numerical order
! by Shell's method (diminishing increment sort)
      integer:: i,j,inc,b
      inc = 1   ! Determine the starting increment
1     inc = 3*inc + 1
      if (inc <= n) GOTO 1
2     continue   ! Loop over the partial sorts
      inc = inc/3
      do i = inc + 1,n   ! Outer loop of straight insertion.
         b = v(i)
         j = i
3        if (v(j - inc) > b) then   ! Inner loop of straight insertion.
            v(j) = v(j - inc)
            j = j - inc
            if (j <= inc) GOTO 4
            GOTO 3
         end if
4        v(j) = b
      end do
      if (inc > 1) GOTO 2
      RETURN
   END SUBROUTINE

   SUBROUTINE sort3(iv)
      integer,intent(inout):: iv(:)
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
      if (iv(3) < iv(2)) call swap(iv(3),iv(2))
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
CONTAINS
      SUBROUTINE swap(x, y)
         integer,intent(inout):: x,y
         integer:: tmp
         tmp = x
         x = y
         y = tmp
      END SUBROUTINE
   END SUBROUTINE

END MODULE

!!>include 'files_mod.f90'

MODULE files_mod
   public get_free_file_unit, myflush
CONTAINS

   SUBROUTINE get_free_file_unit(iu)
      implicit none
      character(len=*),parameter :: sub_name = 'get_free_file_unit'
      integer,intent(out):: iu
      integer,parameter:: istart = 11
      integer,parameter:: unit_max = 100000
      integer:: ios
      logical:: lopen,lexists
      do iu = istart,unit_max
         inquire(unit=iu, iostat = ios, exist = lexists, opened = lopen)
         if (ios == 0) then
            if (lexists .and. .not.lopen) RETURN
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
      open(unit=iu,file=trim(fn),position = 'append')
   END SUBROUTINE myflush

END MODULE files_mod

!!>include 'constants_mod.f90'

MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars, only: angstrom,pi
   implicit none
   real(wp),parameter:: bondl_oh = 0.945_wp/angstrom
   real(wp),parameter:: angle_oh = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: sigma_Obridge = 2.7_wp/angstrom
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
   real(wp),parameter:: sigma_OH = 3.00_wp/angstrom
   real(wp),parameter:: sigma_O = 2.94_wp/angstrom
   real(wp),parameter:: sigma_Si = 4.20_wp/angstrom
   real(wp),parameter:: sigma_C = 2.70_wp/angstrom
   real(wp),parameter:: sigma_N = 3.296_wp/angstrom
   real(wp),parameter:: bondl_O2 = 1.169_wp/angstrom
   real(wp),parameter:: bondl_N2 = 1.097_wp/angstrom
   real(wp),parameter:: bondl_SiO4 = 1.62_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 0.958_wp/angstrom
   real(wp),parameter:: r_tet = bondl_SiO4 + sigma_OH*0.5_wp
   real(wp),parameter:: sigma(0:3) = (/ sigma_OH,sigma_Si,sigma_OH,0.0_wp /)
   real(wp),parameter:: sigma_2(0:3) = sigma*0.5_wp
END MODULE

!!>include 'atom_types.f90'

MODULE atom_types
   USE precision_mod, only: wp
   implicit none
   integer,parameter:: iOxygen = 0
   integer,parameter:: iSilicon = 1
   integer,parameter:: iOxygenH = 2
   integer,parameter:: iHydrogen = 3
   integer,parameter:: ntyp = 3
   character(2),parameter:: atom_name(0:ntyp) = (/ 'O ','Si','OH','H ' /)
   integer,parameter:: ncmax(0:ntyp) = (/2,4,2,1/)
   integer,allocatable:: atom(:)
CONTAINS

   pure FUNCTION name2atom(c)
      integer:: name2atom
      character(*),intent(in):: c
      integer:: i
      do i = 0,ntyp
         if (trim(c) == trim(atom_name(i))) then
            name2atom = i
            exit
         end if
      end do
      if (i > ntyp) name2atom = -1
   END FUNCTION

END MODULE

!!>include 'seaton_mod.f90'

MODULE seaton_mod
   USE precision_mod
   USE global_vars, only: angstrom, K_ev, qstar
   implicit none
! atom types
   integer,parameter:: iO_Sil = 0
   integer,parameter:: iSi_Sil = 1
   integer,parameter:: iO_OH_sil = 2
   integer,parameter:: iH_OH_sil = 3
   integer,parameter:: iH_H2O = 4
   integer,parameter:: iO_H2O = 5
   integer,parameter:: iO_CO2 = 6
   integer,parameter:: iC_CO2 = 7
   integer,parameter:: iN_N2 = 8
   integer,parameter:: iO_O2 = 9
   integer,parameter:: iN_N2charge = 10
   integer,parameter:: iO_O2charge = 11
   integer,parameter:: ntyplj = 9  ! number of LJ types
   real(wp),parameter:: eps_O_sil   = 185.0_wp*K_ev, sig_O_sil  = 2.708_wp/angstrom,q_O_sil = -0.64025_wp*qstar
   real(wp),parameter:: eps_Si_sil  =  0.0_wp*K_ev, sig_Si_sil = 0.0_wp/angstrom,  q_Si_sil = -2.0_wp*q_O_sil
   real(wp),parameter:: eps_O_OH_sil = 185.0_wp*K_ev, sig_O_OH_sil = 3.0_wp/angstrom, q_O_OH_sil = -0.533_wp*qstar
   real(wp),parameter:: eps_H_OH_sil =  0.0_wp*K_ev, sig_H_OH_sil = 0.0_wp/angstrom, q_H_OH_sil =  0.206_wp*qstar
   real(wp),parameter:: eps_H_H2O   =  0.0_wp*K_ev,  sig_H_H2O = 0.0_wp/angstrom,   q_H_H2O =  0.417_wp*qstar
   real(wp),parameter:: eps_O_H2O   = 76.58_wp*K_ev, sig_O_H2O = 3.1506_wp/angstrom,q_O_H2O = -2.0_wp*q_H_H2O
   real(wp),parameter:: eps_O_CO2   = 82.997_wp*K_ev,sig_O_CO2 = 3.064_wp/angstrom, q_O_CO2 = -0.33225_wp*qstar
   real(wp),parameter:: eps_C_CO2   = 29.999_wp*K_ev,sig_C_CO2 = 2.785_wp/angstrom, q_C_CO2 = -2.0_wp*q_O_CO2
   real(wp),parameter:: eps_N_N2    = 34.897_wp*K_ev, sig_N_N2 = 3.3211_wp/angstrom, q_N_N2 = -0.5475_wp*qstar
   real(wp),parameter:: eps_O_O2    = 43.183_wp*K_ev, sig_O_O2 = 3.1062_wp/angstrom, q_O_O2 = -0.3577_wp*qstar
!  bond lengths
   real(wp),parameter:: bondl_SiO = 1.60_wp/angstrom
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
   real(wp),parameter:: bondl_H20 = 0.9572_wp/angstrom
! LJ parameters
   real(wp),parameter:: epsi(0:ntyplj) = (/ &
                        eps_O_sil, &
                        eps_Si_sil, &
                        eps_O_OH_sil, &
                        eps_H_OH_sil, &
                        eps_O_H2O, &
                        eps_H_H2O, &
                        eps_O_CO2, &
                        eps_C_CO2, &
                        eps_N_N2, &
                        eps_O_O2 /)
   real(wp),parameter:: sigi(0:ntyplj) = (/ &
                        sig_O_sil, &
                        sig_Si_sil, &
                        sig_O_OH_sil, &
                        sig_H_OH_sil, &
                        sig_O_H2O, &
                        sig_H_H2O, &
                        sig_O_CO2, &
                        sig_C_CO2, &
                        sig_N_N2, &
                        sig_O_O2 /)
   real(wp),parameter:: sigi2(0:ntyplj) = sigi*0.5_wp
! Charges
   real(wp),parameter:: q_seaton(0:11) = (/ &
                        q_O_sil, &
                        q_Si_sil, &
                        q_O_OH_sil, &
                        q_H_OH_sil, &
                        q_O_H2O, &
                        q_H_H2O, &
                        q_O_CO2, &
                        q_C_CO2, &
                        q_N_N2, &
                        q_O_O2, &
                -2.0_wp*q_N_N2, &
                -2.0_wp*q_O_O2 /)
END MODULE seaton_mod

!!>include 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:)
   real(wp):: boxl,boxli,boxl2
   interface pbc
      module procedure  pbc_v, pbc_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      r(1:2) = r(1:2) - boxl*anint(r(1:2)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      r(:,1:2) = r(:,1:2) - boxl*anint(r(:,1:2)*boxli)
   END SUBROUTINE

   SUBROUTINE WRITE_XYZ(iu,nattached)
      USE atom_types, only: atom_name,atom
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
      read(iu,*) ctmp,ctmp,ctmp,it
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
   USE global_vars, only: natom,nseed
   USE atom_types
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
      if (it /= nattached) stop 'it /= nattached'
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
            RETURN
         end if
         do j = 1,ncmax(atom(K))
            if (proximity(K,j) == M) then
               nearest_neighbor2 = .true.
               RETURN
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
         if (atom(k) == iOxygenH) noh = noh + 1
      end do
   END FUNCTION


   SUBROUTINE check_proximity(consistent,ifirst)
!     Check that if atom i is bonded to k then
!     atom k is bonded to i, if everything is ok
!     then consistent is .true.
      logical,intent(out):: consistent
      integer,intent(out):: ifirst
      integer:: i,j,k
      consistent = .true.
      ifirst = 0
      do i = 1,natom
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (k == 0) cycle
            if ( count(proximity(k,:) == i) /= 1) then
               consistent = .false.
               ifirst = i
               RETURN
            end if
         end do
      end do
   END SUBROUTINE

   SUBROUTINE check_proximity2(consistent,ifirst)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      logical,intent(out):: consistent
      integer,intent(out):: ifirst
      integer:: i,j,k,ia,a1,a2,i1,i2
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
               if (atom(k) == iSilicon) GOTO 100   ! but not Si
               if (atom(k) == iHydrogen) GOTO 100  ! and not H
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = proximity(i,1)
            if (i1 == 0) GOTO 100
            a1 = atom(proximity(i,1))
            if (a1 /= iSilicon) GOTO 100
            if (i > nseed) then  ! non seed bridging O
               i2 = proximity(i,2)  !  must be bonded to  2 Silicons
               if (i2 == 0) GOTO 100
               a2 = atom(proximity(i,2))
               if (a2 /= iSilicon) GOTO 100
            end if
            if (any(proximity(i,3:) /= 0)) GOTO 100
         case(iOxygenH)  ! OH Oxygen must be bonded to Si and H
            i1 = proximity(i,1)
            if (i1 == 0) GOTO 100
            a1 = atom(proximity(i,1))
            if (i > nseed) then
               i2 = proximity(i,2)
               if (i2 == 0) GOTO 100
               a2 = atom(proximity(i,2))
               if (i2 == 0) GOTO 100
               if (.NOT.((a1 == iHydrogen .and. a2 == iSilicon).or. &
                         (a2 == iHydrogen .and. a1 == iSilicon))) GOTO 100
            else
               if (a1 /= iHydrogen) GOTO 100
            end if
            if (any(proximity(i,3:) /= 0)) GOTO 100
         case default
            stop 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ifirst = i
      RETURN
   END SUBROUTINE

   SUBROUTINE delete_atom(i)
!     Delete atom i from the system and renumber the
!     last atom (natom) to number i
      USE coordinates_mod, only: rxyz
      integer,intent(in):: i
      integer:: j,k
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         where( proximity(k,:) == i ) proximity(k,:) = 0
      end do
      rxyz(i,1:3) = rxyz(natom,1:3)
      proximity(i,:) = proximity(natom,:)
      atom(i) = atom(natom)
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
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
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1, -1
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
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1, -1
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
            ni = ni + 1
            ilst(ni) = i
         end do
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1, -1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE

END MODULE connectivity_mod
!!>include 'charges_mod.f90'

MODULE charges_mod
   USE precision_mod
   USE global_vars, only: nseed
   USE atom_types
   USE seaton_mod, only: qi => q_seaton
   implicit none
   PRIVATE
   real(wp),allocatable:: charge(:)
   PUBLIC:: charge, assign_charge, qsi
!
CONTAINS

   pure FUNCTION qsi(n)
      real(wp):: qsi
      integer,intent(in):: n
      qsi = (n - 4)*qi(iOxygen)*0.5_wp - n*(qi(iOxygenH) + qi(iHydrogen))
   END FUNCTION

   SUBROUTINE assign_charge(ifirst,ilast)
      USE connectivity_mod, only: noh
      USE atom_types
      integer,intent(in):: ifirst,ilast
      real(wp):: qo,qh,qoh
      integer:: i
      qoh = qi(iOxygenH)
      qh = qi(iHydrogen)
      qo = qi(iOxygen)
      do i = ifirst,nseed
         select case(atom(i))
         case(iOxygen)
            charge(i) = qo
            charge(i - nseed) = -qo*0.5_wp
         case(iOxygenH)
            charge(i) = qoh
            charge(i - nseed) = -(qoh + qh)
         case default
            stop 'assign_charge: wrong seed type'
         end select
      end do
      do i = nseed + 1,ilast
         select case(atom(i))
         case(iSilicon)
            charge(i) = qsi(noh(i))
         case(iOxygen)
            charge(i) = qo
         case(iOxygenH)
            charge(i) = qoh
         case(iHydrogen)
            charge(i) = qh
         case default
            stop 'assign_charge: unknown atom type'
         end select
      end do
   END SUBROUTINE

END MODULE charges_mod

!!>include 'hydrogens_mod.f90'

MODULE hydrogens_mod
   USE precision_mod
   USE global_vars, only: natom,nseed
   USE constants_mod, only: bondl_oh
   USE atom_types, only: atom,iHydrogen
   USE connectivity_mod, only: proximity
   USE coordinates_mod
   implicit none
CONTAINS
   SUBROUTINE adjust_hydrogens()
      real(wp):: rij(3)
      integer:: j,i
      do j = nseed + 1,natom
         if (atom(j) /= iHydrogen) cycle
         i = proximity(j,1)
         rij(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
         call pbc(rij)
         rij = rij/sqrt(dot_product(rij,rij))
         rxyz(j,1:3) = rxyz(i,1:3) + rij*bondl_oh
         call pbc(rxyz(j,:))
         if (rxyz(j,3) < 0.0_wp) rxyz(j,3) = 0.0_wp
      end do
   END SUBROUTINE

END MODULE hydrogens_mod

!!>include 'bond_list_mod.f90'

MODULE bond_list_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer:: nbondtot
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types, only: atom,ncmax
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

END MODULE bond_list_mod

!!>include 'frames_mod.f90'

MODULE FRAMES_MOD
   USE precision_mod, only: wp
   implicit none
   public new_frame_file, write_frame
CONTAINS

   SUBROUTINE new_frame_file(iu,base_name,iframe)
      USE files_mod
      integer,intent(in):: iframe
      integer,intent(out):: iu
      character(len=*),intent(in):: base_name
      character(len=32):: ctmp
      character(len=132):: frame_file
      integer:: ios
      logical:: lexists,lopen
      call get_free_file_unit(iu)
      write(unit=ctmp,fmt='(i4.4)') iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.pdb'
      inquire(file=trim(frame_file),iostat=ios,exist = lexists,opened = lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit=iu,file=trim(adjustl(frame_file)),form = 'formatted')
   END SUBROUTINE new_frame_file

   SUBROUTINE write_frame(iu,ifirst,ilast)
      USE files_mod
      USE bond_list_mod
      USE global_vars, only: natom,angstrom
      USE connectivity_mod
      USE atom_types, only: atom,atom_name,ncmax
      USE coordinates_mod, only: rxyz,boxl2
      integer,intent(in):: iu,ifirst,ilast
      integer:: i,j,ib,ia,iu2
      logical:: lopen,fix_bonds
      character(len=132):: frame_file
      do i = ifirst,ilast
         write(iu,110)'HETATM',i,atom_name(atom(i)),'   ',i,rxyz(i,1:3)*angstrom
      end do
      call bond_list
!
!     write(iu,*) 'nbondtot = ',nbondtot,' natom = ',natom
!     do i = 1,natom
!        write(iu,*) i,atom_name(atom(i)),nbond(i)
!     end do
!
      inquire(unit=iu,name=frame_file)
      j = len_trim(frame_file)
      frame_file = frame_file(1:j - 3)//'scr'
      fix_bonds = .false.
      do i = 1,nbondtot
         ia = ibond(1,i)
         ib = ibond(2,i)
         if (ANY( ABS(rxyz(ia,1:3) - rxyz(ib,1:3)) > 1.0_wp )) then
            fix_bonds = .true.
            inquire(file=frame_file,opened = lopen)
            if (.not.lopen) then
               call get_free_file_unit(iu2)
               open(unit=iu2,file=trim(frame_file))
            end if
            write(iu2,*)'select ',ia,',',ib
            write(iu2,*)'wireframe off'
         end if
!         write(iu,"(a6,5i5)")'CONECT',ibond(1,i),ibond(2,i)
      end do
      if (fix_bonds) then
         write(iu2,*)
         close(iu2,status='keep')
      end if
      do i = 1,ilast
         if ( all(proximity(i,:) == 0) ) cycle
         write(iu,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,110)'HETATM',99991,'AR','   ',99991, -boxl2*angstrom, -boxl2*angstrom,0.0
      write(iu,110)'HETATM',99992,'AR','   ',99992, -boxl2*angstrom, boxl2*angstrom,0.0
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2*angstrom, -boxl2*angstrom,0.0
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2*angstrom, boxl2*angstrom,0.0
      write(iu,110)'HETATM',99995,'AR','   ',99995, -boxl2*angstrom, -boxl2*angstrom,5.0
      write(iu,110)'HETATM',99996,'AR','   ',99996, -boxl2*angstrom, boxl2*angstrom,5.0
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2*angstrom, -boxl2*angstrom,5.0
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2*angstrom, boxl2*angstrom,5.0
      write(iu,"(a6,5i5)")'CONECT',99991,99992
      write(iu,"(a6,5i5)")'CONECT',99991,99993
      write(iu,"(a6,5i5)")'CONECT',99991,99995
      write(iu,"(a6,5i5)")'CONECT',99992,99994
      write(iu,"(a6,5i5)")'CONECT',99992,99996
      write(iu,"(a6,5i5)")'CONECT',99993,99994
      write(iu,"(a6,5i5)")'CONECT',99993,99997
      write(iu,"(a6,5i5)")'CONECT',99994,99998
      write(iu,"(a6,5i5)")'CONECT',99995,99996
      write(iu,"(a6,5i5)")'CONECT',99995,99997
      write(iu,"(a6,5i5)")'CONECT',99996,99998
      write(iu,"(a6,5i5)")'CONECT',99997,99998
      write(iu,'(a/)')'END'
      close(iu,status='keep')
110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3)
   END SUBROUTINE write_frame

END MODULE FRAMES_MOD

!!>include 'rand_mod.f90'

MODULE rand_mod
! a random number generator (Numerical Recipies)
! with auxillary subroutines for storing and setting the
! state space
   USE precision_mod, only: wp,i4b
   implicit none
   private
   public:: rand,get_rand_state,set_rand_state,read_rand_state,write_rand_state
   integer(i4b),private,save:: ix = -1, iy = -1
   integer(i4b),public,save:: irand_seed
!
CONTAINS

   SUBROUTINE set_rand_state(ix0,iy0)
      integer(i4b),intent(in):: ix0,iy0
      real(wp):: tmp
      tmp = rand()
      ix = ix0
      iy = iy0
   END SUBROUTINE set_rand_state

   SUBROUTINE get_rand_state(ix0,iy0)
      integer(i4b),intent(out):: ix0,iy0
      ix0 = ix
      iy0 = iy
   END SUBROUTINE get_rand_state

   SUBROUTINE read_rand_state(iu)
      integer(i4b),intent(in):: iu
      real(wp):: tmp
      tmp = rand()
      read(iu,*) ix,iy
      rewind(iu)
   END SUBROUTINE read_rand_state

   SUBROUTINE write_rand_state(iu)
      integer(i4b),intent(in):: iu
      write(iu,*) ix,iy
      call flush(iu)
   END SUBROUTINE write_rand_state

   FUNCTION rand()
!-----random number generator
      real(wp):: rand
      integer(i4b),parameter::ia = 16807,im = 2147483647,iq = 127773,ir = 2836
      integer(i4b),save:: k
      real(wp),save::am
      if (irand_seed <= 0 .or. iy < 0 ) then
         am = nearest(1.0_wp, -1.0_wp)/im
         iy = ior(ieor(888889999,abs(irand_seed)),1)
         ix = ieor(777755555,abs(irand_seed))
         irand_seed = abs(irand_seed) + 1
      end if
      ix = ieor(ix,ishft(ix,13))
      ix = ieor(ix,ishft(ix, -17))
      ix = ieor(ix,ishft(ix,5))
      k = iy/iq
      iy = ia*(iy - k*iq) - ir*k
      if (iy < 0) iy = iy + im
      rand = am*ior(iand(im,ieor(ix,iy)),1)
   END FUNCTION rand

! 'minimal standard' RNG
!   FUNCTION rand0()
!      real(wp):: rand0
!      real(wp),parameter:: im=2147483647.0_wp,ia=16807.0_wp
!      real(wp),parameter:: am=1.0_wp/im
!      irand_seed = mod(ia*mod(ia*irand_seed,im),im)
!      rand0 = am*irand_seed
!   END FUNCTION

END MODULE rand_mod

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
      crossp_3d(1) = v1(2)*v2(3) - v1(3)*v2(2)
      crossp_3d(2) = v1(3)*v2(1) - v1(1)*v2(3)
      crossp_3d(3) = v1(1)*v2(2) - v1(2)*v2(1)
   END FUNCTION crossp_3d

   PURE FUNCTION transpose_3d(m)
      real(wp) :: transpose_3d(3,3)
      real(wp),intent(in) :: m(3,3)
      integer :: i,j
      forall(i = 1:3,j = 1:3) transpose_3d(j,i) = m(i,j)
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

!!>include 'rotate_axis_mod.f90'

MODULE rotate_axis_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE zaxis2vect(r1,aa)
!     Return the rotation matrix aa that rotates the
!     z-axis to the unit vector r1
      real(wp),intent(in):: r1(3)
      real(wp),intent(out):: aa(3,3)
      real(wp),parameter:: eps = epsilon(0.0_wp)
      real(wp):: a,b,c,a11,a12,a22,a2pb2
      a = r1(1)
      b = r1(2)
      c = r1(3)
      a2pb2 = a*a + b*b
      if ( a2pb2 > eps ) then
         a11 = (b*b + a*a*c)/a2pb2
         a22 = (1.0_wp + c) - a11   ! (a*a + b*b*c)/(a*a + b*b)
         a12 = -a*b/(1.0_wp + c)    ! (a*b*(c - 1))/(a*a + b*b)
         aa(1,1:3) = (/ a11, a12, a /)
         aa(2,1:3) = (/ a12, a22, b /)
         aa(3,1:3) = (/  -a,  -b, c /)
      else
         aa = 0.0_wp
         if ( c > 0.0_wp ) then
            aa(1,1) = 1.0_wp
            aa(2,2) = 1.0_wp
            aa(3,3) = 1.0_wp
         else
            aa(2,1) = -1.0_wp
            aa(1,2) = -1.0_wp
            aa(3,3) = -1.0_wp
         end if
      end if
   END SUBROUTINE zaxis2vect

END MODULE rotate_axis_mod

!!>include 'insert_seeds_mod_H.f90'

MODULE insert_seeds_mod
   USE precision_mod
   implicit none
CONTAINS

   SUBROUTINE insert_seeds()
      USE rand_mod
      USE global_vars, only: natom,pi
      USE constants_mod, only: bondl_oh,angle_oh,bondl_SiO4,sigma_OH
      USE coordinates_mod
      USE connectivity_mod
      USE atom_types
      integer,parameter:: ntrymax = 1000000
      real(wp):: xx,yy,zz,r(3),cangle,sangle,phi
      integer:: itry,j,i
      cangle = cos(angle_oh)
      sangle = sin(angle_oh)
      natom = 0
! Try to place nseed Oxygen (OH) atoms on the z=0 surface
! ensuring no overlaps
      seed_loop: do itry = 1,ntrymax
         xx = boxl*rand() - boxl2
         yy = boxl*rand() - boxl2
         zz = 0.0_wp
         do j = 1,natom
            r(1) = xx - rxyz(j,1)
            r(2) = yy - rxyz(j,2)
            r(3) = zz - rxyz(j,3)
            call pbc(r)
            if (dot_product(r,r) < (sigma_OH)**2) CYCLE seed_loop
         end do
         natom = natom + 1
         rxyz(natom,1:3) = (/ xx,yy,zz /)
         atom(natom) = iOxygenH
         if (natom == nseed) EXIT
      end do seed_loop
      if (itry >= ntrymax) then
         nseed = natom
         write(*,*) 'actual no. seeds used = ',nseed
      end if
! place 'dummy' Si atoms below surface for electrostatic neutrality
      forall(i = 1:nseed) rxyz(i - nseed,:) = rxyz(i,:) - (/0.0_wp,0.0_wp,bondl_SiO4/)
      atom(1 - nseed:0) = iSilicon
! attach a Hydrogen to each Oxygen
! overlaps allowed: H positions will be optimised later
      do i = 1,nseed
         natom = nseed + i
         proximity(i,1) = natom
         proximity(natom,1) = i
         atom(natom) = iHydrogen
         phi = rand()*2.0_wp*pi
         rxyz(natom,1) = rxyz(i,1) + bondl_oh*sangle*cos(phi)
         rxyz(natom,2) = rxyz(i,2) + bondl_oh*sangle*sin(phi)
         rxyz(natom,3) = rxyz(i,3) - bondl_oh*cangle
         call pbc(rxyz(natom,:))
      end do
   END SUBROUTINE

END MODULE insert_seeds_mod

!!>include 'nlist_mod.f90'

MODULE NLIST_MOD
      USE precision_mod, only: wp
      USE global_vars, only: natom_max,natom
      USE coordinates_mod, only: rxyz,boxl2
      implicit none
      PRIVATE
      integer,public:: neigh
      real(wp),private:: delx,dely,delz,delmin
      real(wp),private:: delxi,delyi,delzi
      integer:: ncelx,ncely,ncelz,ncelt
      integer,parameter:: ncelmax = 100000,neighmx = 125
      integer,public:: hoc(ncelmax),hoc_old(ncelmax),ncell(neighmx)
      integer,allocatable,public:: ll_old(:),ll(:),lr(:)
      public:: NEIGCELL,cell,INIT_NLIST,NEW_NLIST,print_cell
      public:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll
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

   pure FUNCTION nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   END FUNCTION

   SUBROUTINE INIT_NLIST(Lx,Ly,Lz,Rc)
      real(wp),intent(in):: Lx,Ly,Lz,Rc
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max))
      ncelx = INT(Lx/Rc)
      ncely = INT(Ly/Rc)
      ncelz = INT(Lz/Rc)
      ncelt = ncelx*ncely*ncelz
      delx = Lx/ncelx
      dely = Ly/ncely
      delz = Lz/ncelz
      delmin = min(delx,dely,delz)
      delxi = 1.0_wp/delx
      delyi = 1.0_wp/dely
      delzi = 1.0_wp/delz
      !write(*,*)'ncelx = ',ncelx
      !write(*,*)'ncely = ',ncely
      !write(*,*)'ncelz = ',ncelz
      !write(*,*)'ncelt = ',ncelt
      !write(*,*)'delx = ',delx
      !write(*,*)'dely = ',dely
      !write(*,*)'delz = ',delz
      LL(:) = 0
      HOC(:) = 0
      LR(:) = 0
   END SUBROUTINE

   PURE FUNCTION CELL(XYZ)
!     determines cell number for position x,y,z
      real(wp),intent(in)::  XYZ(:)
      integer:: CELL
      CELL = INT((XYZ(1) + boxl2)*delxi) &
           + INT((XYZ(2) + boxl2)*delyi)*ncelx &
           + INT( XYZ(3)       *delzi)*ncelx*ncely + 1
      RETURN
   END FUNCTION CELL

   SUBROUTINE NEW_NLIST
!     makes a new neighbour list using the linked-list algorithm
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head - of - chain
      ! make linked list:
      do i = 1,natom
         ! determine cell number
         ic = CELL(rxyz(i,:))
         ! update linked - list and head of chain
         LL(i) = HOC(ic)
         if (HOC(ic) /= 0) LR(HOC(ic)) = i
         HOC(ic) = i
      end do
      RETURN
   END SUBROUTINE NEW_NLIST

   SUBROUTINE PUSH_CELL(i,ic)
      integer,intent(in):: i,ic
      LL(i) = HOC(ic)
      if (HOC(ic) /= 0) LR(HOC(ic)) = i
      HOC(ic) = i
   END SUBROUTINE

   SUBROUTINE POP_CELL(i,ic)
      integer,intent(in):: i,ic
      if (HOC(ic) == i) then
         HOC(ic) = LL(i)
         if (LL(i) /= 0) LR(LL(i)) = 0
         LL(i) = 0
      else
         LL(LR(i)) = LL(i)
         if (LL(i) /= 0) LR(LL(i)) = LR(i)
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
      if (mod(ic,n2) == 0) icz = icz - 1
      icy = (ic - (icz - 1)*n2)/ncelx + 1
      if (mod(ic - (icz - 1)*n2,ncelx) == 0) icy = icy - 1
      icx = ic - (icy - 1)*ncelx - (icz - 1)*n2
      ncell = 0
      do iz = -nlayer,nlayer
         iccz = icz + iz
         if (iccz < 1) then
            CYCLE ! NOTE: NOT Periodic in Z - direction
         else if (iccz > ncelz) then
            CYCLE ! NOTE: NOT Periodic in Z - direction
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
               inn = iccx + (iccy - 1)*ncelx + (iccz - 1)*n2
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
      !icy = (ic - (icz - 1)*n2)/ncelx + 1
      !icx = ic - (icy - 1)*ncelx - (icz - 1)*n2
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
      if (mod(ic,(ncelx*ncely)) == 0) icz = icz - 1
      icy = (ic - (icz - 1)*(ncelx*ncely))/ncelx + 1
      if (mod(ic - (icz - 1)*(ncelx*ncely),ncelx) == 0) icy = icy - 1
      icx = ic - (icy - 1)*ncelx - (icz - 1)*(ncelx*ncely)
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
!            if (nc == 0) cycle
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
      forall(i = 1:n2) ncell(i) = n2*(iz - 1) + i
   END SUBROUTINE

END MODULE NLIST_MOD

!!>include 'verlet_list_nl.f90'

MODULE VERLET_LIST
   USE precision_mod, only: wp
   USE global_vars, only: natom_max,natom
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nabormx = 50
   real(wp),allocatable:: rxyz_old(:,:)
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
      call NEW_NLIST
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

END MODULE VERLET_LIST

!!>include 'tetra_coords_mod.f90'
MODULE tetra_coords_mod
      USE precision_mod
      implicit none
CONTAINS

   SUBROUTINE tetra_coords(bondl,tetra)
      real(wp),intent(in):: bondl
      real(wp),intent(out):: tetra(3,4)
!
! returns a set of coordinates for the 4 outer atoms of a
! tetrahedral molecule (with the central atom at the origin).
! Atom 1 is on the z axis below the origin and atoms 2,3, & 4
! are on a plane above the origin.
!
      real(wp):: zz,sth,cth,rxy,xx,yy
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,2) = (/ xx, yy, zz /)
      tetra(1:3,3) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,4) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

   SUBROUTINE tetra_coords5(bondl,tetra)
      real(wp),intent(in):: bondl
      real(wp),intent(out):: tetra(3,5)
!
! returns a set of coordinates for a tetrahedral molecule
! with the central atom at the origin (atom 1), one atom
! on the z axis below it (atom 2), and atoms 3,4, and 5 on
! a plane above the origin.
!
      real(wp):: zz,sth,cth,rxy,xx,yy
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

   SUBROUTINE SiOH4_coords(bondl,bondl_OH,angle_OH,tetra)
      USE rotate_axis_mod
      USE rand_mod
      real(wp),intent(in):: bondl,bondl_OH,angle_OH
      real(wp),intent(out):: tetra(3,9)
!
! Returns a set of coordinates for a tetrahedral molecule
! with the central atom at the origin (atom 1), one atom
! on the z axis below it (atom 2), and atoms 3,4, and 5 on
! a plane above the origin.
! Atoms 6,7,8 and 9 are attached to atoms 2,3,4 and 5
! respectively
!
      real(wp):: zz,sth,cth,rxy,xx,yy,rr(3)
      real(wp):: aa(3,3),cang,sang,pi,phi
      integer:: j
      pi = 4.0_wp*(atan(1.0_wp))
      cang = cos(angle_oh)
      sang = sin(angle_oh)
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
! attach atoms 6-9 to 2-5
      do j = 2,5
         rr = tetra(1:3,j)/bondl
         call zaxis2vect(rr,aa)
         phi = rand()*2.0_wp*pi
         zz = -bondl_oh*cang
         xx = bondl_oh*sang*cos(phi)
         yy = bondl_oh*sang*sin(phi)
         tetra(1:3,j + 4) = tetra(1:3,j) + matmul(aa,(/ xx,yy,zz /))
      end do
   END SUBROUTINE

END MODULE

!!>include 'quat2mat_mod.f90'
MODULE quat2mat_mod
   USE precision_mod
   implicit none
CONTAINS

   PURE SUBROUTINE quat2mat(Q,M)
      real(wp),intent(in):: Q(4)
      real(wp),intent(out):: M(3,3)
      real(wp):: tx,ty,tz,tq,tk
! q is a quaternion describing a random rotation
! m is the rotation matrix
      TX = Q(2)*Q(2)
      TY = Q(3)*Q(3)
      TZ = Q(4)*Q(4)
      TQ = TY + TZ
      if ((TQ + TX + Q(1)*Q(1)) /= 0.0_wp) then
         TK = 2.0_wp / (TQ + TX + Q(1)*Q(1))
      else
         TK = 0.0_wp
      end if
      M(1,1) = 1.0_wp - TK*TQ
      M(2,2) = 1.0_wp - TK*(TX + TZ)
      M(3,3) = 1.0_wp - TK*(TX + TY)
      TX = TK*Q(2)
      TY = TK*Q(3)
      TQ = (TK*Q(4))*Q(1)
      TK = TX*Q(3)
      M(1,2) = TK - TQ
      M(2,1) = TK + TQ
      TQ = TY*Q(1)
      TK = TX*Q(4)
      M(1,3) = TK + TQ
      M(3,1) = TK - TQ
      TQ = TX*Q(1)
      TK = TY*Q(4)
      M(2,3) = TK - TQ
      M(3,2) = TK + TQ
   END SUBROUTINE quat2mat

END MODULE quat2mat_mod

!!>include 'ran_point_sphere_mod.f90'

MODULE ran_point_sphere
   USE precision_mod
   USE rand_mod, only: rand
   implicit none
CONTAINS

   SUBROUTINE ran_point_in_sphere(qx,qy,qz)
      real(wp),intent(out):: qx,qy,qz
      real(wp)::s1
      do
         qx = (2.0_wp*rand() - 1.0_wp)
         qy = (2.0_wp*rand() - 1.0_wp)
         qz = (2.0_wp*rand() - 1.0_wp)
         s1 = qx*qx + qy*qy + qz*qz
         if (s1 <= 1.0_wp) exit
      end do
      RETURN
   END SUBROUTINE

   SUBROUTINE ran3sph(x,y,z)
!
! random vector on the surface of a sphere
! Marsaglia,G., Ann.math.Stat.,43,645-6,(1972)
!
      real(wp),intent(out):: x,y,z
      real(wp):: s,tmp

1     continue
         x = 2.0_wp*rand() - 1.0_wp
         y = 2.0_wp*rand() - 1.0_wp
         s = x*x + y*y
      if (s >= 1.0_wp) GOTO 1
      z = 2.0_wp*s - 1.0_wp
      tmp = 2.0_wp*sqrt(1.0_wp - s)
      x = x*tmp
      y = y*tmp
      RETURN
   END SUBROUTINE

   SUBROUTINE ran3sph2(x,y,z)
      real(wp),parameter::twopi = 6.283185307179586476925286766559005768394_wp
      real(wp),intent(out)::x,y,z
      real(wp)::s,tmp
      z = 2.0_wp*rand() - 1.0_wp
      tmp = rand()*twopi
      s = sqrt(1.0_wp - z*z)
      x = s*cos(tmp)
      y = s*sin(tmp)
      RETURN
   END SUBROUTINE

   SUBROUTINE ran_point_on_sphere(qx,qy,qz)
      real(wp),intent(out):: qx,qy,qz
      real(wp)::s1,r1
      do
        qx = (2.0_wp*rand() - 1.0_wp)
        qy = (2.0_wp*rand() - 1.0_wp)
        qz = (2.0_wp*rand() - 1.0_wp)
        s1 = qx*qx + qy*qy + qz*qz
        if (s1 <= 1.0_wp) exit
      end do
      r1 = 1.0_wp/sqrt(s1)
      qx = qx*r1
      qx = qx*r1
      qz = qx*r1
      RETURN
   END SUBROUTINE

   SUBROUTINE ran4sph(q)
! random point on surface of 4-Sphere (a Quaternion)
      real(wp),intent(out)::q(4)
      real(wp)::s1,s2
      do
        q(1) = (2.0_wp*rand() - 1.0_wp)
        q(2) = (2.0_wp*rand() - 1.0_wp)
        s1 = q(1)*q(1) + q(2)*q(2)
        if (s1 <= 1.0_wp) exit
      end do
      do
        q(3) = (2.0_wp*rand() - 1.0_wp)
        q(4) = (2.0_wp*rand() - 1.0_wp)
        s2 = q(3)*q(3) + q(4)*q(4)
        if (s2 <= 1.0_wp) exit
      end do
      q(3) = q(3)*sqrt((1.0_wp - s1)/s2)
      q(4) = q(4)*sqrt((1.0_wp - s1)/s2)
      RETURN
   END SUBROUTINE

END MODULE

!include 'Keating_parameters_Burlakov.f90'
!include 'repul_energy_VL_Burlakov.f90'
!!>include 'Keating_parameters_vonAlfthan.f90'

MODULE Keating_parameters
!
! Parameters for the simplified Keating potential.
! Source: S. von alfthan et al., Phys. Rev. B 68, 073203, (2003)
!
   USE precision_mod, only: wp
   implicit none
   real(wp),parameter:: KSiSi = 908.0_wp   ! eV nm^-2
   real(wp),parameter:: KSiO = 2700.0_wp
   real(wp),parameter:: ASiSi = 0.235_wp   ! nm
   real(wp),parameter:: ASiO = 0.161_wp
   real(wp),parameter:: KSiSiSi = 3.58_wp  ! eV
   real(wp),parameter:: KSiSiO = 3.93263_wp
   real(wp),parameter:: KOSiO = 4.32_wp
   real(wp),parameter:: KSiOSi = 2.0_wp
   real(wp),parameter:: kbond(1:3) = (/ KSiO,KSiSi,KSiO /)
   real(wp),parameter:: abond(1:3) = (/ ASiO,ASiSi,ASiO /)
   real(wp),parameter:: acosang(0:2) = (/ 1.0_wp,1.0_wp/3.0_wp,1.0_wp /)
END MODULE

!!>include 'repul_energy_VL_vonAlfthan.f90'

MODULE repul_energy_mod
   USE precision_mod, only: wp
   implicit none
   PUBLIC :: repul_energy
CONTAINS

   PURE FUNCTION repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE VERLET_LIST
      real(wp):: repul_energy
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2
      integer:: M,L,jj
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repul_energy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
         if (atom(L) == iHydrogen) CYCLE outer_atom_loop
         atom_loop: do jj = 1, NLIST(L)  ! Verlet - list
            M = LIST(L,jj)
            if (atom(M) == iHydrogen) CYCLE atom_loop
            if (nearest_neighbor2(L,M)) CYCLE atom_loop
            if (M == L) CYCLE atom_loop
            rr(1:3) = rxyz(M,1:3) - rxyz(L,1:3)
            call pbc(rr)
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            if (r2 <= rden2) repul_energy = repul_energy + 0.5_wp*ff*(r2 - rden2)**2
         end do atom_loop
      end do outer_atom_loop
   END FUNCTION repul_energy

END MODULE repul_energy_mod

!!>include 'Keating.f90'

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
      kang(0,1,0) = KOSiO
      kang(2,1,0) = KOSiO
      kang(0,1,2) = KOSiO
      kang(2,1,2) = KOSiO
      !! kang(0,1,1) = KSiSiO  ! Not used for now (no Si - Si bonds)
      !! kang(1,1,0) = KSiSiO
      !! kang(1,1,1) = KSiSiSi
      kang(1,0,1) = KSiOSi
      kang(1,2,1) = KSiOSi
      ebond = 0.0_wp
      eangle = 0.0_wp
!
      do i = ifirst,ilast
      if (atom(i) == iHydrogen) cycle
!     bond energy
      do L = 1,ncmax(atom(i))
         j = proximity(i,L)
         if (atom(j) == iHydrogen) cycle
         if (j /= 0) then
         rij(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
         call pbc(rij)
         ebond = ebond + 0.5_wp*kbond(atom(i) + atom(j)) &
               *(sqrt(dot_product(rij,rij)) - abond(atom(i) + atom(j)))**2
!!             *(len_3d(rij) - abond(atom(i)+atom(j)))**2
         end if
      end do
!
!     angle energy
      do L = 1,ncmax(atom(i)) - 1
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

!!>include 'HKNonLattice.f90'

MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster = 0,n_cluster_old = 0
   integer,parameter,private :: ncmx = 100
   integer,allocatable:: atomL(:),atomL_old(:)
   integer,allocatable,private:: cluster(:,:)
   integer,private:: clusterC(ncmx)
CONTAINS

   SUBROUTINE Init_HKNonLattice(natom_max)
      integer,intent(in):: natom_max
      allocate(atomL(natom_max),atomL_old(natom_max),cluster(ncmx,natom_max))
      atomL = 0
      atomL_old = 0
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
      integer:: i,ii,j,nnmax,iNodeNextL,inn,NodeLPmin,k
      integer:: RelabL1(NumberOfNodes),RelabL(NumberOfNodes),RelabLB(NumberOfNodes)
      !
      ! STEPS 1,2 & 3 of AL - Futaisi and Tadeusz Patzek algorithm
      !
      nnmax = size(NodeNext,dim=2)
      NumberOfClusters = 0
      NodeL(:) = 0
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
               NodeNextL(iNodeNextL) = NodeL(j)
            end do
            if ( any(NodeNextL > 0) ) then  ! Case 4c ii: a labeled neighbor exists
               ! Put in the minimum labeling
               inn = 0
               do ii = 1,iNodeNextL
                  j = NodeNextL(ii)
                  if (j == 0) cycle
                  inn = inn + 1
                  N(inn) = NodeLP(j)
               end do
               NodeLPmin = minval(NodeLP(N(1:inn)))
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
      where(NodeLP1(2:iNodeLP) > NodeLP1(1:iNodeLP - 1)) RelabLB(2:iNodeLP) = 1
      RelabL1(1:iNodeLP - 1) = NodeLP1(2:iNodeLP)*RelabLB(2:iNodeLP)
      RelabL = 0
      RelabL(1) = NodeLP1(1)
      ii = 1
      do i = 1,iNodeLP - 1
         if (RelabL1(i) == 0) cycle
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
      RETURN
   END SUBROUTINE HKNonLattice

! some routines used mainly for debugging
   SUBROUTINE analyse_cluster(natom)
      implicit none
      integer,intent(in):: natom
      integer:: i,L
      clusterC(1:n_cluster) = 0
      cluster(1:n_cluster,1:ncmx) = 0
      do i = 1,natom
         L = AtomL(i)
         clusterC(L) = clusterC(L) + 1
         cluster(L,clusterC(L)) = i
      end do
      ! some consistency checks
      do i = 1,n_cluster
         if (clusterC(i) /= count(atomL(1:natom) == i)) then
            stop 'error in analyse_cluster'
         end if
      end do
      if (sum(clusterC(1:n_cluster)) /= natom) then
         stop 'error in analyse_cluster: wrong count'
      end if
   END SUBROUTINE analyse_cluster

   SUBROUTINE print_cluster(natom,iu)
      implicit none
      integer,intent(in):: natom,iu
      integer:: i
      do i = 1,n_cluster
         write(iu,*)'cluster = ',i,'count = ',clusterC(i)
         write(iu,'(10i6)') cluster(i,1:clusterC(i))
         if (clusterC(i) /= count(atomL(1:natom) == i)) then
            stop 'error in print_cluster'
         end if
      end do
   END SUBROUTINE print_cluster

END MODULE HKNonLattice_mod

!!>include 'probe_mol_mod.f90'

MODULE probe_mol_mod
   USE precision_mod
   USE constants_mod
   implicit none
!
!  definition of probe_mol derived type
!  the allocatable component requires a standard extension to F95
!  described in the technical report ISO/IEC TR 15581
!
   TYPE probe_mol
      integer:: n
      real(wp),allocatable:: r(:,:)
      integer,allocatable:: atom(:)
      real(wp),allocatable:: rad(:),q(:)
   END TYPE probe_mol
!
   type(probe_mol):: probe_H2O, probe_SiO4, probe_SiOH4
   type(probe_mol):: probe_N2, probe_O2, probe_tet, probe_CO2
!
CONTAINS

   SUBROUTINE new_probe_mol(P,n,radius,atomtype)
      type(probe_mol),intent(out):: P
      integer,intent(in):: n
      real(wp),intent(in),optional:: radius
      integer,intent(in),optional:: atomtype
      P%n = n
      allocate(P%r(3,n),P%atom(n),P%rad(n),P%q(n))
      P%r = 0.0_wp
      if (present(radius)) P%rad = radius
      if (present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types
      USE seaton_mod
      USE charges_mod
      call new_probe_mol(probe_H2O,1,sigma_H2O*0.5_wp)
      call new_probe_mol(probe_tet,1,r_tet)
      call new_probe_mol(probe_O2,2,sigma_O*0.5_wp)
      call new_probe_mol(probe_N2,2,sigma_N*0.5_wp)
      call new_probe_mol(probe_CO2,3,sigma_O*0.5_wp)
      call new_probe_mol(probe_SiO4,4,sigma_OH*0.5_wp)
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO4,bondl_OH,angle_OH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigi2(probe_SiOH4%atom)
      probe_SiOH4%q = q_seaton(probe_SiOH4%atom)
      probe_SiOH4%q(1) = qsi(4)
      call tetra_coords(bondl_SiO4,probe_SiO4%r)
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigi2(probe_SiO4%atom)
      probe_SiO4%q = q_seaton(probe_SiO4%atom)

      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_O2*0.5_wp /)

      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_N2*0.5_wp /)
      probe_CO2%rad(1) = sigma_C*0.5_wp
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp, -bondl_CO2 /)
      write(*,*) 'SiOH4'
      write(*,*) probe_SiOH4%n
      write(*,'(9f12.6)') probe_SiOH4%rad
      write(*,'(9f12.6)') probe_SiOH4%q
      write(*,'(3f12.6)') probe_SiOH4%r
      write(*,*) 'SiO4'
      write(*,*) probe_SiO4%n
      write(*,'(9f12.6)') probe_SiO4%rad
      write(*,'(3f12.6)') probe_SiO4%r
      write(*,*) 'O2'
      write(*,*) probe_O2%n
      write(*,'(9f12.6)') probe_O2%rad
      write(*,'(3f12.6)') probe_O2%r
      write(*,*) 'N2'
      write(*,*) probe_N2%n
      write(*,'(9f12.6)') probe_N2%rad
      write(*,'(3f12.6)') probe_N2%r
      write(*,*) 'CO2'
      write(*,*) probe_CO2%n
      write(*,'(9f12.6)') probe_CO2%rad
      write(*,'(3f12.6)') probe_CO2%r
   END SUBROUTINE init_probe_mols

   SUBROUTINE WRITE_PROBE_MOL(iu,P)
      USE atom_types
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      integer:: i
      do i = 1,P%n
         write(iu,110)'HETATM',i,atom_name(P%atom(i)),'   ',i,P%r(1:3,i)*angstrom
      end do
110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3)
   END SUBROUTINE

END MODULE probe_mol_mod

!!>include 'voidage4_mod.f90'

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
      integer:: k,j,ic,jj,nc,it,nat
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
         probe(1:2,:) = probe(1:2,:) - anint(probe(1:2,:)*boxli)*boxl
! check for overlap
         probe_atom_loop: do k = 1,nat
            ic = CELL( probe(1:3,k) )
            call NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj = 1,neigh
               nc = ncell(jj)
               if (nc == 0) cycle cell_loop
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
            if (nc == 0) cycle cell_loop
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

!!>include 'deposition4_mod_H.f90'

MODULE deposition_mod
   USE precision_mod, only: wp
   implicit none
   logical,parameter:: use_voidage_calc = .TRUE.
   integer,parameter,private:: ntrial = 10000
   integer,parameter,private:: nbinmax = 1000
   integer:: nbin,itop_nonp,itop_nonp_SiO4,itop_nonp_H2O
   real(wp):: nden(nbinmax,0:3),voidage(nbinmax)
   real(wp):: intvoidage(nbinmax),delz
CONTAINS

   elemental FUNCTION voidage_approx(no,sigmap)
      real(wp):: voidage_approx
      real(wp),intent(in):: no,sigmap
      real(wp),parameter:: alpha = 34.7_wp
      real(wp),parameter:: ro = 1.5_wp
      real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
      real(wp):: t1,t2
      t1 = -(4.0_wp/3.0_wp)*pi*no*(ro*ro*ro)
      t2 = 1.0_wp + alpha*no*(sigmap/(2.0_wp*ro))
      voidage_approx = exp(t1*t2*t2*t2)
   END FUNCTION voidage_approx

   SUBROUTINE den_profile(ibegin,iend,dl)
      USE coordinates_mod
      USE atom_types
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,delzi
      integer:: j,ii,jj
      top_atom = maxval(rxyz(ibegin:iend,3))
      nbin = int(top_atom/dl) + 1
      delz = dl
      delzi = 1.0_wp/delz
      nden(1:nbin,:) = 0.0_wp
      do j = ibegin,iend
         ii = int(rxyz(j,3)*delzi + 1.0_wp)
         if (ii > nbin) then
            if ( rxyz(j,3) == top_atom) then
               ii = nbin
            else
               stop 'ii > nbin'
            end if
         end if
         jj = atom(j)
         nden(ii,jj) = nden(ii,jj) + 1.0_wp
      end do
      nden(1:nbin,:) = nden(1:nbin,:)/(boxl*boxl*delz*angstrom**3)
   END SUBROUTINE den_profile

   SUBROUTINE integrate_voidage(ibegin,iend,fv,ifv)
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: fv(:)
      real(wp),intent(out):: ifv(:)
      real(wp):: tot
      integer:: i
      tot = 0.0_wp
      do i = ibegin,iend
         tot = tot + fv(i)
         ifv(i) = tot
      end do
      ifv(ibegin:iend) = ifv(ibegin:iend)/tot
   END SUBROUTINE

   SUBROUTINE voidage_profile(probe,dl)
      USE voidage_mod
      USE coordinates_mod
      USE probe_mol_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: dl
      real(wp):: top_atom,zl,zu
      integer:: i
      top_atom = maxval(rxyz(1:natom,3))
      nbin = int(top_atom/dl) + 1
      delz = dl
      do i = 1,nbin
         zl = (i - 1)*delz
         zu = i*delz
         if (probe%n == 1) then
            voidage(i) = voidage_calc1(ntrial,zl,zu,probe%rad(1))
         else
            voidage(i) = voidage_calc(ntrial,zl,zu,probe)
         end if
      end do
   END SUBROUTINE

   SUBROUTINE prob_dist(void_crit)
      real(wp),intent(in):: void_crit
      integer:: i
      itop_nonp = 1
      do i = nbin,1, -1
         if ( voidage(i) < void_crit ) then
            itop_nonp = i + 1
            exit
         end if
      end do
      call integrate_voidage(itop_nonp,nbin,voidage,intvoidage)
   END SUBROUTINE

   SUBROUTINE write_den_profile(iu)
      integer,intent(in):: iu
      integer:: i
      write(iu,'(a)') '#nden1 atoms/[A^3]'
      write(iu,'(a,i4)') '#itop_nonp = ',itop_nonp
      do i = 1,nbin
         write(iu,'(5f12.6)') (i - 0.5_wp)*delz,nden(i,:)
      end do
      write(iu,'(/)')
   END SUBROUTINE

   SUBROUTINE write_voidage_dist(iu)
      integer,intent(in):: iu
      integer:: i
      real(wp):: rsum
      rsum = sum(voidage(itop_nonp:nbin))
      write(iu,'(a)') '# voidage'
      do i = itop_nonp,nbin
         write(iu,'(3f12.6)') (i - 0.5_wp)*delz,voidage(i)/rsum,intvoidage(i)
      end do
      write(iu,'(/)')
   END SUBROUTINE

   SUBROUTINE write_voidage(iu)
      integer,intent(in):: iu
      integer:: i
      write(iu,'(a)') '# voidage'
      do i = 1,nbin
         write(iu,'(2f12.6)') (i - 0.5_wp)*delz,voidage(i)
      end do
      write(iu,'(/)')
      call flush(iu)
   END SUBROUTINE

   SUBROUTINE rand_pos(z_dep)
      USE rand_mod
      real(wp),intent(out):: z_dep
      integer:: iv1(1)
      real(wp):: rnd
      rnd = rand()
      iv1 = MINLOC(ABS(intvoidage(itop_nonp:nbin) - rnd)) + itop_nonp - 1
!      z_dep = (iv1(1)-0.5)*delz
      z_dep = (iv1(1) - 1)*delz + rand()*delz
   END SUBROUTINE

   SUBROUTINE sample_z_bin(iz)
      USE rand_mod
      integer,intent(out):: iz
      integer:: iv1(1)
      real(wp):: rnd
      rnd = rand()
      iv1 = MINLOC(ABS(intvoidage(itop_nonp:nbin) - rnd)) + itop_nonp - 1
      iz = iv1(1)
   END SUBROUTINE

   SUBROUTINE sample_z_bin_uniform(iz)
      USE rand_mod
      integer,intent(out):: iz
      real(wp):: rnd
      rnd = rand()
      iz = int(rnd*(nbin - itop_nonp + 1)) + itop_nonp
   END SUBROUTINE

END MODULE deposition_mod

!!>include 'rdf_mod.f90'

MODULE RDF_MOD
   USE precision_mod, only: wp
   USE global_vars, only: natom,nseed
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nbin = 200
   real(wp),private:: RDF(nbin)
   real(wp),private:: DL1
   PUBLIC:: RDF_CALC, RDF_PRINT
CONTAINS

   SUBROUTINE RDF_CALC
      real(wp):: rv(3),r,r2,rmax,rmax2
      integer:: i,j,k
      rmax = boxl*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDF = 0.0_wp
      do i = nseed + 1,natom - 1
         do j = i + 1,natom
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rmax2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1) + 1
            RDF(k) = RDF(k) + 1.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC

   SUBROUTINE RDF_PRINT(io,nattached)
      integer,intent(in):: io,nattached
      integer:: i
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom - nseed) = ',natom - nseed
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i - 0.5_wp)*DL1,RDF(i)/(natom - nseed)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT

END MODULE RDF_MOD

!!>include 'relax2_VL_mod.f90'

MODULE relax_mod
   USE precision_mod, only: wp
   implicit none
   public:: RELAX
CONTAINS

   SUBROUTINE RELAX(nr,ifirst,ilast,delt,alist,na,quench,scale_repul)
      USE coordinates_mod
      USE atom_types
      USE connectivity_mod
      USE Keating
      USE repul_energy_mod
      USE verlet_list
      USE rand_mod
      USE global_vars, only: etot
      implicit none
      integer,intent(in):: nr,ifirst,ilast,alist(:),na
      logical,intent(in):: quench,scale_repul
      real(wp),intent(in):: delt
      real(wp),parameter:: de = 0.006_wp
!     real(wp):: de = 0.08_wp
      integer:: ir,i,j,k,nfail
      real(wp):: coprxyz(3),energy1,energy2,del_energy
      real(wp):: EEE,fscale0,fscale,E0
      logical:: special_atom,any_special_atoms
!
      nfail = 0
      any_special_atoms = ( alist(1) /= 0 )
      do ir = 1,NR

      if (scale_repul) fscale0 = real(ir - 1,wp)/real(NR - 1,wp)
      if (quench) E0 = Etot + (NR - ir)*delt/(NR - 1)

      do i = ifirst,ilast
         if (atom(i) == iHydrogen) CYCLE
         if (any_special_atoms) then
            special_atom = any(alist(1:na) == i)
         else
            special_atom = .false.
         end if
         if (special_atom .and. scale_repul) then
            fscale = fscale0
         else
            fscale = 1.0_wp
         end if

         if (update(i)) call new_vlist
         coprxyz(1:3) = rxyz(i,1:3)
         energy1 = energy4(i,i)
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (atom(k) == iHydrogen) CYCLE
            if (k > 0) then
               energy1 = energy1 + energy4(k,k)
            end if
         end do
         energy1 = energy1 + fscale*repul_energy(i,i)

         do
            rxyz(i,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
            if (rxyz(i,3) > 0.0_wp) EXIT
         end do
         rxyz(i,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(i,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
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
         energy2 = energy2 + fscale*repul_energy(i,i)
!
         del_energy = energy2 - energy1

         if (special_atom) then
            if (quench) then
               EEE = del_energy/E0
            else
               EEE = del_energy/delt
            end if
         else
            EEE = del_energy/Etot
         end if

         if (EEE <  0.0_wp) GOTO 22
         if (EEE > 50.0_wp) GOTO 11
         if (rand() < exp(-EEE)) GOTO 22
11       CONTINUE
            rxyz(i,1:3) = coprxyz(1:3)
            nfail = nfail + 1
22       CONTINUE
      end do

      end do
!     print '(f0.6,2(i7,1x),f0.6)',de,NR*(ilast-ifirst+1),nfail,real(nfail)/(NR*(ilast-ifirst+1))
!     if (real(nfail)/(NR*(ilast-ifirst+1)) < 0.5) then
!        de = de*1.01
!     else
!        de = de/1.01
!     end if
   END SUBROUTINE RELAX

END MODULE relax_mod

!!>include 'attachment_NL_mod_H.f90'

MODULE attachment_mod
   USE precision_mod
   USE global_vars
   USE constants_mod
   USE atom_types
   USE coordinates_mod
   USE connectivity_mod
   USE rand_mod
   implicit none
CONTAINS

   SUBROUTINE get_atom_for_attachment4(iatom)
      integer,intent(out):: iatom
      integer:: arraybonds(natom_max),n,L
      n = 0
      arraybonds(1:natom) = 0
      do L = 1,natom
         if (atom(L) /= iOxygenH) cycle
         n = n + 1
         arraybonds(n) = L
      end do
      iatom = arraybonds( int(rand()*n) + 1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment3(iz,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      integer,intent(out):: iatom
      integer:: arraybonds(natom_max),n,L,ic,nc
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0) cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            n = n + 1
            arraybonds(n) = L
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      iatom = arraybonds( int(rand()*n) + 1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment2(iz,diam,bondl,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      real(wp),intent(in):: diam,bondl
      integer,intent(out):: iatom
      integer,parameter:: ntrial_max = 100
      real(wp):: r3v(3),alfa,beta
      integer:: arraybonds(natom_max),i,n,L,it,ic,nc
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0) cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            trial_loop: do it = 1,ntrial_max
               alfa = acos(2.0_wp*rand() - 1.0_wp)
               beta = rand()*2.0_wp*pi
               rxyz(natom + 1,3) = rxyz(L,3) + bondl*cos(beta)
               rxyz(natom + 1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
               rxyz(natom + 1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
               call pbc(rxyz(natom + 1,:))
               do i = 1,natom
                  if (atom(i) == iHydrogen) cycle
                  if (atom(i) == iSilicon) cycle
                  if (i == L) cycle
                  r3v(1:3) = rxyz(i,1:3) - rxyz(natom + 1,1:3)
                  call pbc(r3v)
                  if (dot_product(r3v,r3v) <= diam**2) CYCLE trial_loop
               end do
               exit trial_loop
            end do trial_loop
            if (it < ntrial_max) then
               n = n + 1
               arraybonds(n) = L
            end if
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      iatom = arraybonds( int(rand()*n) + 1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment(z_dep,diam,bondl,iatom)
      real(wp),intent(in):: z_dep,diam,bondl
      integer,intent(out):: iatom
      integer,parameter:: ntrial_max = 100
      real(wp):: r3v(3),alfa,beta
      integer:: arraybonds(natom_max),i,n,L,it,iv(1)
!
      n = 0
      arraybonds(1:natom) = 0
      atom_loop: do L = 1,natom
         if (atom(L) /= iOxygenH) cycle atom_loop
         trial_loop: do it = 1,ntrial_max
            alfa = acos(2.0_wp*rand() - 1.0_wp)
            beta = rand()*2.0_wp*pi
            rxyz(natom + 1,3) = rxyz(L,3) + bondl*cos(beta)
            rxyz(natom + 1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
            rxyz(natom + 1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
            call pbc(rxyz(natom + 1,:))
            do i = 1,natom
               if (atom(i) == iHydrogen) cycle
               if (atom(i) == iSilicon) cycle
               if (i == L) cycle
               r3v(1:3) = rxyz(i,1:3) - rxyz(natom + 1,1:3)
               call pbc(r3v)
               if (dot_product(r3v,r3v) <= diam**2) CYCLE trial_loop
            end do
            exit trial_loop
         end do trial_loop
         if (it >= ntrial_max) cycle atom_loop
         n = n + 1
         arraybonds(n) = L
      end do atom_loop
      iv = MINLOC(ABS(rxyz(arraybonds(1:n),3) - z_dep))
      iatom = arraybonds(iv(1))
      if (iatom <= nseed) then
         do
            iatom = rand()*nseed + 1.0
            if (atom(iatom) == iOxygenH) exit
         end do
      end if
   END SUBROUTINE get_atom_for_attachment

   SUBROUTINE attach_tetrahedra(L,beta,bondl,success)
      USE matvec3d_mod,only: len_3d,matvec_3d
      USE nlist_mod
      USE HKNonLattice_mod
      USE rotate_axis_mod
      integer,intent(in):: L
      real(wp),intent(in):: beta,bondl
      logical,intent(out):: success
      integer,parameter:: ntry = 10
      real(wp):: xx,yy,zz,phi,cphi,sphi,cbeta,sbeta
      real(wp):: sth,cth,rxy,r2,cangle,sangle
      real(wp):: r3v(3),r12(3),r23(3),r34(3),r35(3),r36(3),rh(3)
      real(wp):: tetra(3,3),aa(3,3),aa1(3,3),aa2(3,3),aa3(3,3)
      integer:: i,j,LSi,ii,jj,LS,M,ic,nc,nl,LH,itry,ilist(7)
      logical:: overlap
      nl = nlayers_ll(0.3_wp)
      cbeta = cos(beta)
      sbeta = sin(beta)
      cangle = cos(angle_oh)
      sangle = sin(angle_oh)
      success = .false.
!-----Attaching a new Silicon atom
      if (L > nseed) then
         do j = 1,ncmax(atom(L))
            ii = proximity(L,j)
            if (ii == 0) then
               STOP 'attach_tetrahedra: ERROR 1'
            end if
            if ( atom(ii) == iSilicon ) then
               LSi = ii
            else if (atom(ii) == iHydrogen) then
               LH = ii
            else
               STOP 'attach_tetrahedra: ERROR 2'
            end if
         end do
         r3v(1:3) = rxyz(L,1:3) - rxyz(LSi,1:3)
         call pbc(r3v)
         r12 = r3v/len_3d(r3v)
      else ! attach to a seed
         LH = proximity(L,1)
         r12 = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      end if
      rh = rxyz(LH,:)  ! store the Hydrogen position
      LS = LH          ! the new Silicon replaces the Hydrogen
      call zaxis2vect( r12, aa )
      do itry = 1,ntry
         phi = rand()*2.0_wp*pi
         cphi = cos(phi)
         sphi = sin(phi)
         zz = -bondl*cbeta
         xx = bondl*sbeta*cphi
         yy = bondl*sbeta*sphi
         r23 = matmul(aa,(/xx,yy,zz/))
         rxyz(LS,1:3) = rxyz(L,1:3) + r23
         if (rxyz(LS,3) > 0.05_wp) exit ! Silicon position is ok
         ! need to try again
         if (itry == ntry) then ! give up
            rxyz(LH,:) = rh
            RETURN
         end if
      end do
!
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      r23 = r23/len_3d(r23)
      call zaxis2vect(r23,aa)
!
!     attach the 3 Oxygen atoms
      do itry = 1,ntry
         phi = rand()*2.0_wp*pi
         cphi = cos(phi)
         sphi = sin(phi)
         xx = rxy*sphi
         yy = -rxy*cphi
         tetra(1:3,1) = (/ xx, yy, zz /)
         tetra(1:3,2) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
         tetra(1:3,3) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
!
         r34 = matmul(aa,tetra(1:3,1))
         r35 = matmul(aa,tetra(1:3,2))
         r36 = matmul(aa,tetra(1:3,3))
         rxyz(natom + 1,1:3) = rxyz(LS,1:3) + r34
         rxyz(natom + 2,1:3) = rxyz(LS,1:3) + r35
         rxyz(natom + 3,1:3) = rxyz(LS,1:3) + r36
! if the new atom is below the surface (z = 0) try again
         if (ANY(rxyz(natom + 1:natom + 3,3) <= 0.0_wp)) then
            if (itry == ntry) then
               rxyz(LH,:) = rh
               RETURN
            else
               cycle
            end if
         else
            exit
         end if
      end do
!
! attach Hydrogens to Oxygens
      r34 = r34/len_3d(r34)
      r35 = r35/len_3d(r35)
      r36 = r36/len_3d(r36)
      call zaxis2vect(r34,aa1)
      call zaxis2vect(r35,aa2)
      call zaxis2vect(r36,aa3)
      zz = -bondl_oh*cangle
      do itry = 1,ntry
         phi = rand()*2.0_wp*pi
         xx = bondl_oh*sangle*cos(phi)
         yy = bondl_oh*sangle*sin(phi)
         rxyz(natom + 4,1:3) = rxyz(natom + 1,1:3) + matmul(aa1,(/ xx,yy,zz /))
         phi = rand()*2.0_wp*pi
         xx = bondl_oh*sangle*cos(phi)
         yy = bondl_oh*sangle*sin(phi)
         rxyz(natom + 5,1:3) = rxyz(natom + 2,1:3) + matmul(aa2,(/ xx,yy,zz /))
         phi = rand()*2.0_wp*pi
         xx = bondl_oh*sangle*cos(phi)
         yy = bondl_oh*sangle*sin(phi)
         rxyz(natom + 6,1:3) = rxyz(natom + 3,1:3) + matmul(aa3,(/ xx,yy,zz /))
! if the new atom is below the surface (z = 0) try again
         if (ANY(rxyz(natom + 4:natom + 6,3) <= 0.0_wp)) then
            if (itry == ntry) then
               rxyz(LH,:) = rh
               RETURN
            else
               cycle
            end if
         else
            exit
         end if
      end do
!
!-----Apply PBC to 7 added atoms
      natom = natom + 6
      call pbc(rxyz(LS,:))
      call pbc(rxyz(natom - 5:natom,:))
!
      call set_proximity(L,LH,LS)
      atom(L) = iOxygen
      atom(LS) = iSilicon
      proximity(LS,1) = L
      proximity(LS,2) = natom - 5
      proximity(LS,3) = natom - 4
      proximity(LS,4) = natom - 3
      proximity(natom - 5,1) = LS
      proximity(natom - 4,1) = LS
      proximity(natom - 3,1) = LS
      proximity(natom - 5,2) = natom - 2
      proximity(natom - 4,2) = natom - 1
      proximity(natom - 3,2) = natom
      proximity(natom - 2,1) = natom - 5
      proximity(natom - 1,1) = natom - 4
      proximity(natom  ,1) = natom - 3
      atom(natom - 5:natom - 3) = iOxygenH
      atom(natom - 2:natom) = iHydrogen
!
      atomL(natom - 5:natom) = atomL(L)  ! cluster label for new atoms
      call NEW_NLIST
!
! Check for overlap
!!    ilist = (/ natom-5,natom-4,natom-3,natom-2,natom-1,natom /)
      ilist(1:3) = (/ natom - 5,natom - 4,natom - 3 /)
      overlap = .false.
      outer_atom_loop: do i = 1,3
      ii = ilist(i)
      ic = CELL(rxyz(ii,:))  ! link list
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if (M == ii) GOTO 100
            if (M > natom - 3) GOTO 100  ! ignore the last 3 Hydrogens
            if (nearest_neighbor2(ii,m)) GOTO 100
            r3v(1:3) = rxyz(M,1:3) - rxyz(ii,1:3)
            call pbc(r3v)
            r2 = r3v(1)**2 + r3v(2)**2 + r3v(3)**2
            if (r2 <= (sigma_2(atom(m)) + sigma_2(atom(ii)))**2) then
               overlap = .true.
               exit outer_atom_loop
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
!
      if (overlap) then
         rxyz(LH,:) = rh
         atom(LH) = iHydrogen
         proximity(LH,:) = 0
         proximity(LH,1) = L
         atom(L) = iOxygenH
         atom(natom - 5:natom) = 0
         proximity(natom - 5:natom,:) = 0
         natom = natom - 6
         call NEW_NLIST
         RETURN
      else
         success = .true.
         RETURN
      end if
   END SUBROUTINE

END MODULE attachment_mod

!!>include 'add_water2_VL_mod_H.f90'

MODULE add_water_mod
   USE precision_mod, only: wp
   implicit none
   real(wp):: fwater
   real(wp),save:: sumfwater = 0.0_wp
   integer,save:: iadd_water = 0
   integer,save:: ntry_add_water
CONTAINS

   SUBROUTINE GET_NUM_ADD_WATER_ATTEMPTS()
      if (fwater >= 1.0_wp) then
         ntry_add_water = int(fwater)
         RETURN
      else
         sumfwater = sumfwater + fwater
         if (sumfwater >= 1.0_wp) then
            sumfwater = 0.0_wp
            ntry_add_water = 1
         else
            ntry_add_water = 0
         end if
      end if
   END SUBROUTINE

   SUBROUTINE add_water(iz,success)
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE constants_mod, only: pi,angle_oh,bondl_oh
      USE Relax_mod
      USE Keating
      USE repul_energy_mod
      USE verlet_list
      USE nlist_mod
      USE rand_mod
      USE HKNonLattice_mod
      USE ran_point_sphere,only: ran3sph
      USE deposition_mod
      USE matvec3d_mod,only: len_3d,matvec_3d
      USE rotate_axis_mod
      integer,intent(in):: iz
      logical,intent(out):: success
      integer,parameter:: ntry = 10
      integer:: arraybonds(natom_max),copatom(natom_max)
      integer:: copproximity(natom_max,size(proximity,dim=2))
      real(wp):: coprxyz(natom_max,3)
      real(wp):: EEE,energy1,energy2 !!,sx,sy,sz,rad
      integer:: iat,i,noxyb,iv(1),j,k,ir,itmp
      integer:: natom_old,iSi1,iSi2,iSi,ns,itry
      real(wp):: Eden,del_energy,z_pos,xx,yy,zz,phi,cangle,sangle
      real(wp):: aa1(3,3),aa2(3,3),r1(3),r2(3)
      logical :: possible,dangling_group !!,hot_atoms
      success = .FALSE.
      z_pos = (iz - 1)*delz + rand()*delz
      cangle = cos(angle_oh)
      sangle = sin(angle_oh)
!
! find bridging Oxygens
! do not consider seed oxygens (for now)
!
      noxyb = 0
      do i = 2*nseed + 1,natom
         if (atom(i) /= iOxygen) CYCLE
         noxyb = noxyb + 1
         arraybonds(noxyb) = i
      end do
!
! Now consider only the seed oxygen atoms which are bonded to
! a single Si(OH)3 group
!
      seed_loop: do i = 1,nseed
         if (atom(i) /= iOxygen) cycle
         ! Pick iSi , the Silicon attached to i
         iSi = 0
         do j = 1,ncmax(atom(i))
            iSi = proximity(i,j)
            if (iSi /= 0) exit
         end do
         if (.NOT.is_dangling_group(iSi,i)) cycle seed_loop
         noxyb = noxyb + 1
         arraybonds(noxyb) = i
      end do seed_loop

      if (noxyb == 0) RETURN

      Eden = Etot
!!    hot_atoms=.false.
!
!     pick the bridging oxygen closest to z_pos
3131  CONTINUE
      possible = .false.
      dangling_group = .FALSE.
      iv = MINLOC(ABS( rxyz(arraybonds(1:noxyb),3) - z_pos ))
      iat = arraybonds(iv(1))
!     but, if it's a seed pick one of the suitable seeds at random
      if (iat <= nseed) then
         ns = count(arraybonds(1:noxyb) <= nseed)
         if (ns > 1) then
            ir = rand()*nseed + 1
            k = 0
            do i = 1,noxyb
               if (arraybonds(i) <= nseed) k = k + 1
               if (k >= ir) then
                  iat = arraybonds(i)
                  exit
               end if
            end do
         end if
      end if
!
      if (iat > nseed) then
         ! The two attached Silicons
         iSi1 = proximity(iat,1)
         iSi2 = proximity(iat,2)
! if there is a dangling Si(OH)3 group the Silicon atom is iSi1
         if (is_dangling_group(iSi1,iat)) then
            dangling_group = .TRUE.
         end if
         if (is_dangling_group(iSi2,iat)) then
            if (dangling_group) STOP 'add_water: ERROR, both groups cannot be dangling'
            dangling_group = .TRUE.
            itmp = iSi2
            iSi2 = iSi1
            iSi1 = itmp
         end if
      else
         dangling_group = .TRUE.
         iSi1 = maxval(proximity(iat,:))
      end if
      if (dangling_group) then
         if (iz < itop_nonp_SiO4) then
            arraybonds(iv) = arraybonds(noxyb)
            arraybonds(noxyb) = 0
            noxyb = noxyb - 1
            if (noxyb > 0) GOTO 3131
            write(*,*)'add_water: fail - water attack not allowed'
            RETURN
         end if
      end if
!
! save all the old data before the attempted move
      call NEW_VLIST
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
      natom_old = natom
      do i = 1,natom
         coprxyz(i,1:3) = rxyz(i,1:3)
         copproximity(i,:) = proximity(i,:)
         copatom(i) = atom(i)
         cnlist(i) = nlist(i)
         clist(i,1:nlist(i)) = list(i,1:nlist(i))
         atomL_old(i) = atomL(i)
      end do
!
! energy before the move
      energy1 = energy4(1,natom) + repul_energy(1,natom)
! add the water
      atom(iat) = iOxygenH
      natom = natom + 1
! now: natom = O, natom+1 = H1, natom+2 = H2
      atom(natom) = iOxygenH
      atom(natom + 1) = iHydrogen
      atom(natom + 2) = iHydrogen
      atomL(natom) = atomL(iSi1)  ! cluster label
      atomL(natom + 1) = atomL(iSi1)
      atomL(natom + 2) = atomL(iat)
      call set_proximity(natom,  0, iSi1)
      call set_proximity(natom,  0, natom + 1)
      call set_proximity(iSi1, iat, natom)
      call set_proximity(iat, iSi1, natom + 2)
      proximity(natom + 1,1) = natom
      proximity(natom + 2,1) = iat
      natom = natom + 2
! now: natom = H2, natom-1 = H1, natom-2 = O
!
      if (dangling_group) then
         possible = .TRUE.
      else
         ! First check if structure is still connected
         call HKNonLattice(natom,proximity,n_cluster,atomL)
         if (n_cluster > nseed) GOTO 1020
         if ( .NOT.ANY( atomL(1:nseed) == atomL(iSi1) ) ) GOTO 1020
         if ( .NOT.ANY( atomL(1:nseed) == atomL(iSi2) ) ) GOTO 1020
         possible = .true.
      end if
!     positions for the new atoms
      rxyz(natom - 2,:) = rxyz(iat,:)
      r1(1:3) = rxyz(natom - 2,1:3) - rxyz(iSi1,1:3)
      call pbc(r1)
      r1 = r1/len_3d(r1)
      if (iat > nseed) then
         r2(1:3) = rxyz(iat,1:3) - rxyz(iSi2,1:3)
         call pbc(r2)
         r2 = r2/len_3d(r2)
      else
         r2 = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      end if
      call zaxis2vect(r1,aa1)
      call zaxis2vect(r2,aa2)
      zz = -bondl_oh*cangle
! trying to add the 2 Hydrogen atoms
      do itry = 1,ntry
         phi = rand()*2.0_wp*pi
         xx = bondl_oh*sangle*cos(phi)
         yy = bondl_oh*sangle*sin(phi)
         rxyz(natom - 1,1:3) = rxyz(natom - 2,1:3) + matmul(aa1,(/ xx,yy,zz /))
         phi = rand()*2.0_wp*pi
         xx = bondl_oh*sangle*cos(phi)
         yy = bondl_oh*sangle*sin(phi)
         rxyz(natom,1:3) = rxyz(iat,1:3) + matmul(aa2,(/ xx,yy,zz /))
! if the new atom is below the surface (z = 0) try again
         if (ANY(rxyz(natom - 1:natom,3) <= 0.0_wp)) then
            if (itry == ntry) then
               write(*,*) 'add_water: could not place Hydrogens above surface'
               GOTO 1020
            else
               cycle
            end if
         else
            exit
         end if
      end do
! apply PBC
      call pbc(rxyz(natom - 1,:))
      call pbc(rxyz(natom,:))
!
! move the Oxygens (iat and natom) to reasonable (?) starting positions
!
!!      do
!!         call ran3sph(sx,sy,sz)
!!         rad = sigma_2(atom(iat))
!!         rxyz(natom,1:3) = rxyz(natom,1:3) + (/sx,sy,sz/)*rad
!!         if (rxyz(natom,3) > 0.0_wp) exit
!!      end do
!!      if (iat > nseed) then
!!         rxyz(iat,1:3) = rxyz(iat,1;3) - (/sx,sy,sz/)*rad
!!         if (rxyz(iat,3)<=0.0_wp) rxyz(iat,3) =0.0_wp
!!      end if
!
      call NEW_VLIST
!
!     call RELAX(NRELAX, nseed+1, natom, hot_atoms)
      call RELAX(NRELAX, nseed + 1, natom, delt = e_activ*0.5_wp, &
                 alist = (/natom,iat/),na = 2, &
                 quench = .TRUE., scale_repul = .TRUE.)
      energy2 = energy4(1,natom) + repul_energy(1,natom)
      del_energy = energy2 - energy1 + del_rxn
      EEE = del_energy/Eden
      if (EEE <  0.0_wp) GOTO 1019
      if (EEE > 50.0_wp) GOTO 1020
      if (rand() < exp(-EEE)) GOTO 1019

1020  CONTINUE ! the move failed so undo the changes
      rxyz(natom_old + 1:natom_old + 3,:) = 0.0_wp
      proximity(natom_old + 1:natom_old + 3,:) = 0
      atomL(natom_old + 1:natom_old + 3) = 0
      natom = natom_old
      n_cluster = n_cluster_old
      do i = 1,natom
         rxyz(i,1:3) = coprxyz(i,1:3)
         proximity(i,:) = copproximity(i,:)
         atom(i) = copatom(i)
         nlist(i) = cnlist(i)
         list(i,1:cnlist(i)) = clist(i,1:cnlist(i))
         atomL(i) = atomL_old(i)
      end do
      call NEW_NLIST
      if (.not.possible) then
         arraybonds(iv) = arraybonds(noxyb)
         arraybonds(noxyb) = 0
         noxyb = noxyb - 1
         if (noxyb > 0) GOTO 3131
         write(*,*) 'add_water: fail - connectivity'
      else
         write(*,*) 'add_water: fail - energy'
      end if
      RETURN

1019  CONTINUE
      if (dangling_group) then
         call delete_group2(iSi1)
         call NEW_VLIST
      end if
      write(*,*) 'add_water: success'
      write(79,*) nattached,rxyz(iat,3)
      success = .TRUE.
      RETURN
   END SUBROUTINE add_water

END MODULE add_water_mod

!!>include 'relax_dangling2_VL_NL_mod_H.f90'

MODULE relax_dangling_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE relax_dangling(nouter)
      USE global_vars
      USE constants_mod, only: bondl_oh
      USE coordinates_mod
      USE connectivity_mod
      USE Relax_mod
      USE Keating
      USE repul_energy_mod
      USE rand_mod
      USE HKNonLattice_mod
      USE verlet_list
      USE nlist_mod
      USE matvec3d_mod,only: len_3d
      USE sort_mod
      integer,intent(in):: nouter
      integer:: arraybonds(natom_max),arraySi(natom_max),sarraySi(natom_max)
      integer:: copproximity(natom_max,size(proximity,dim=2)),copatom(natom_max)
      real(wp):: coprxyz(natom_max,3)
      real(wp):: r3v(3),EEE,energy1,energy2
      integer:: iouter,iat,iat2,i,N_OH_groups,i1,iSi1,ii,nSi,iSi2,ilst(3),ilstO(3)
      integer:: j,L1,s3,jj,k,natom_old,nc,ic,nl,iatH,iat3
      real(wp),parameter:: r = 0.4_wp
      real(wp),parameter:: r2 = r*r
      real(wp):: Eden,del_energy
      logical :: hot_atoms,remove_h2o,switch
!
      nl = nlayers_ll(r)
      if (nouter <= 3) then ! the last 3 added Oxygens are natom - 3, natom - 4, natom - 5
         ilstO(1) = natom - 3
         ilstO(2) = natom - 4
         ilstO(3) = natom - 5
      end if
      outer_loop: do iouter = 1,nouter

!     Pick iat, the Oxygen with a dangling bond
      if (nouter > 3) then
         Eden = Etot
         hot_atoms = .false.
         N_OH_groups = 0
         do i = 1,natom
            if (atom(i) == iOxygenH) then
               N_OH_groups = N_OH_groups + 1
               arraybonds(N_OH_groups) = i
            end if
         end do
         s3 = int(rand()*N_OH_groups) + 1
         iat = arraybonds(s3)
         if (N_OH_groups == 0) EXIT outer_loop
      else
         iat = ilstO(iouter)
         Eden = 0.5_wp
         hot_atoms = .true.
         if (atom(iat) /= iOxygenH) then
            print *,'ERROR: relax_dangling'
            print *,'iat = ',iat,' type = ',atom(iat)
            print *,atom_name(atom(iat))
            STOP
         end if
      end if

!     Pick iSi1, the Silicon attached to iat
      iSi1 = 0
      do i = 1,ncmax(atom(iat))
         ii = proximity(iat,i)
         if (atom(ii) == iSilicon) then
            iSi1 = ii
         else if (atom(ii) == iHydrogen) then
            iatH = ii
         end if
      end do
! if the DB was on a seed
      if (iSi1 == 0) iSi1 = iat
! Look for a Silicon with 'OH' attached and  within distance r of iSi1
      nSi = 0
      ic = CELL(rxyz(iSi1,:))
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
      nc = ncell(jj)
      if (nc == 0) cycle cell_loop
      j = HOC(nc)
      do while (j /= 0)
         if (atom(j) /= iSilicon) GOTO 100
         if (j == iSi1) GOTO 100
         do ii = 1,ncmax(atom(j))
            k = proximity(j,ii)
            if (k == 0) cycle
            if (atom(k) == iOxygenH) exit
            if (ii == ncmax(atom(j))) GOTO 100
         end do
         r3v(1:3) = rxyz(j,1:3) - rxyz(iSi1,1:3)
         call pbc(r3v(:))
         if ((r3v(1)**2 + r3v(2)**2 + r3v(3)**2) <= r2) then
            nSi = nSi + 1
            arraySi(nSi) = j
         end if
100      j = LL(j)
      end do
      end do cell_loop
!
 call qsort(nSi,arraySi(1:nSi),sarraySi(1:nSi))
      if (nSi < 1) CYCLE outer_loop

! save all the old data before the attempted move
      call NEW_VLIST
      do i = 1,natom
         coprxyz(i,1:3) = rxyz(i,1:3)
         copproximity(i,:) = proximity(i,:)
         copatom(i) = atom(i)
         cnlist(i) = nlist(i)
         clist(i,1:nlist(i)) = list(i,1:nlist(i))
         atomL_old(i) = atomL(i)
      end do
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
      natom_old = natom
!
! energy before the move
      energy1 = energy4(1,natom) + repul_energy(1,natom)

! Choose one of the Silicons with a dangling bond
      L1 = int(rand()*nSi) + 1
      iSi2 = sarraySi(L1)
!!    iSi2 = arraySi(L1)
      do
         i1 = int(rand()*ncmax(atom(iSi2))) + 1
         iat2 = proximity(iSi2,i1)
         if (iat2 /= 0) exit
         if (iat2 == 0) then
            write(*,*)'iat2 = ',iat2
            write(*,*)'i1 = ',i1
            write(*,'(a,i6,a,6i5)') 'iSi2 = ',iSi2,': ',proximity(iSi2,1:ncmax(atom(iSi2)))
            write(*,'(a,i6,a,6i5)') 'iSi1 = ',iSi1,': ',proximity(iSi1,1:ncmax(atom(iSi1)))
            write(*,'(a,i6,a,6i5)') 'iat  = ',iat,  ': ',proximity(iat,1:ncmax(atom(iat)))
            stop 'iat2 == 0'
         end if
      end do
      iat3 = 0
      if (iat2 > nseed) then
         do i = 1,ncmax(atom(iat2))
            ii = proximity(iat2,i)
            if (ii == iSi2) cycle
            if (ii == 0) cycle
            iat3 = ii
         end do
      end if
      atom(iat) = iOxygen
      atom(iat2) = iOxygenH
      proximity(iSi2,i1) = iat
      call set_proximity(iat,  iatH, iSi2)
      call set_proximity(iat2, iSi2, iatH)
      call set_proximity(iatH,iat,iat2)
!
      remove_h2o = .false.
      switch = .true.
      if ((iat2 > nseed).AND.atom(iat3) == iHydrogen) then
         switch = .false.
         remove_h2o = .true.
         ilst = (/iatH,iat3,iat2/)
         call shell(3,ilst)
         call delete_atom(ilst(3))
         call delete_atom(ilst(2))
         call delete_atom(ilst(1))
         call new_vlist
      else
         r3v = rxyz(iatH,1:3) - rxyz(iat2,1:3)
         rxyz(iatH,1:3) = rxyz(iat2,1:3) + bondl_oh*r3v/len_3d(r3v)
      end if
!
!     First check if structure is still connected
!
      call HKNonLattice(natom,proximity,n_cluster,atomL)

      if (n_cluster > nseed) GOTO 1020
      if (.NOT.remove_h2o) then
         if ( .NOT.ANY( atomL(1:nseed) == atomL(iat2) ) ) GOTO 1020
      end if
      if ( .NOT.ANY( atomL(1:nseed) == atomL(iSi2) ) ) GOTO 1020
!
!     call RELAX(NRELAX, nseed+1, natom, hot_atoms)
      if (hot_atoms) then
         if (switch) then
            call RELAX(NRELAX, nseed + 1, natom, delt = 0.5_wp, &
                       alist = (/natom - 3,natom - 2,natom - 1,natom/),na = 4, &
                       quench = .TRUE., scale_repul = .FALSE.)
         else
            call RELAX(NRELAX, nseed + 1, natom, delt = e_activ*0.5_wp, &
                       alist = (/natom - 3,natom - 2,natom - 1,natom,iat/),na = 5, &
                       quench = .TRUE., scale_repul = .FALSE.)
         end if
      else
         if (switch) then
            call RELAX(NRELAX, nseed + 1, natom, delt = 0.0_wp, &
                    alist = (/0/),na = 1, &
                    quench = .FALSE., scale_repul = .FALSE.)
         else
            call RELAX(NRELAX, nseed + 1, natom, delt = e_activ*0.5_wp, &
                    alist = (/iat/),na = 1, &
                    quench = .TRUE., scale_repul = .FALSE.)
         end if
      end if
!
      energy2 = energy4(1,natom) + repul_energy(1,natom)
      del_energy = energy2 - energy1
      EEE = del_energy/Eden
      if (EEE <  0.0_wp) GOTO 1019
      if (EEE > 50.0_wp) GOTO 1020
      if (rand() < exp(-EEE)) GOTO 1019
1020  CONTINUE ! reject the changes
      natom = natom_old
      n_cluster = n_cluster_old
      do i = 1,natom
         rxyz(i,1:3) = coprxyz(i,1:3)
         proximity(i,:) = copproximity(i,:)
         atom(i) = copatom(i)
         nlist(i) = cnlist(i)
         list(i,1:cnlist(i)) = clist(i,1:cnlist(i))
         atomL(i) = atomL_old(i)
      end do

      call NEW_NLIST
      cycle outer_loop

1019  CONTINUE
      if (remove_h2o) then
         write(*,*) 'H2O removal ',rxyz(iat,3)
         write(80,*) nattached,rxyz(iat,3)
      else
         write(*,*) 'bond switch ',rxyz(iat,3)
         write(81,*) nattached,rxyz(iat,3)
      end if

      end do outer_loop
      RETURN
   END SUBROUTINE relax_dangling

END MODULE relax_dangling_mod

!!>include 'surface_atoms_mod.f90'

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
         rad = rprobe + sigma_2(atom(ia))
         r0(1:3) = rxyz(ia,1:3)
         trial_loop: do it = 1,ntrial
            call ran3sph(sx,sy,sz)
            r(1:3) = r0(1:3) + (/sx,sy,sz/)*rad
            ic = CELL(r)
            CALL NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj = 1,neigh
               nc = ncell(jj)
               if (nc == 0) cycle cell_loop
               j = HOC(nc)
               do while (j /= 0)
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
      fsurfa = real(count(surface(ibegin:iend)),wp)/(iend - ibegin + 1)
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
         if (nc == 0) cycle cell_loop
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
      write(iu,'(5f12.6)') (iz - 0.5_wp)*delz,fsur,fsur/voidage(iz)
      end do
   END SUBROUTINE

END MODULE SURFACE_ATOMS_MOD

!!>include 'lj_el_mod.f90'

MODULE lj_el_mod
    USE precision_mod, only: wp
    USE global_vars
    implicit none
    private
    save
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: vdwcut(:,:),vdwsf(:,:)
    real(wp),parameter:: rcut = 12.0_wp/angstrom
    real(wp),parameter:: rcut2 = rcut**2
    PUBLIC:: LJ_INIT,ENERGY_LJ_EL
CONTAINS

   SUBROUTINE LJ_INIT
      USE seaton_mod, only: epsi,sigi,ntyplj
      real(wp):: eij,rstij
      integer:: i,j
      if (.not.allocated(aij)) allocate( aij(0:ntyplj,0:ntyplj) )
      if (.not.allocated(bij)) allocate( bij(0:ntyplj,0:ntyplj) )
      if (.not.allocated(vdwcut)) allocate( vdwcut(0:ntyplj,0:ntyplj) )
      if (.not.allocated(vdwsf))  allocate(  vdwsf(0:ntyplj,0:ntyplj) )
      do j = 0,ntyplj    ! lj aij,bij parameters & shifted terms
      do i = 0,ntyplj
        eij = sqrt(epsi(i)*epsi(j))
        rstij = (sigi(i) + sigi(j))*0.5
        aij(i,j) = 4.0*(rstij**12)*eij
        bij(i,j) = 4.0*(rstij**6)*eij
        vdwcut(i,j) = rcut**6*(aij(i,j)*rcut**6 - bij(i,j))
        vdwsf(i,j) = -rcut**6*(12.0*aij(i,j)*rcut**6 - 6.0*bij(i,j))*rcut
        write(*,'(2i6,4g18.9)') i,j,aij(i,j),bij(i,j),vdwcut(i,j),vdwsf(i,j)
      end do
      end do
   END SUBROUTINE

   PURE SUBROUTINE ENERGY_LJ_EL(prxyz,patom,pcharge,nat,sysfirst,syslast,Ulj,Uel)
      USE coordinates_mod
      USE atom_types, only: atom
      USE charges_mod, only: charge
      real(wp),intent(in):: prxyz(:,:),pcharge(:)
      integer,intent(in):: patom(:),nat,sysfirst,syslast
      real(wp),intent(out):: Ulj,Uel
      real(wp):: rr(3),r2,r1,rr6,dele,qL
      integer:: L,M,iaL
      Ulj = 0.0_wp
      Uel = 0.0_wp
      outer_atom_loop: do L = 1,nat
      iaL = patom(L)
      qL = pcharge(L)
      atom_loop: do M = sysfirst,syslast
         rr(1:3) = rxyz(M,1:3) - prxyz(1:3,L)
         rr(1:2) = rr(1:2) - boxl*anint(rr(1:2)*boxli)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         !!call pbc(rr)
         !!r2 = dot_product(rr,rr)
         r1 = sqrt(r2)
         if (r2 < rcut2 .and. M > 0) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(iaL,atom(M))*rr6 - bij(iaL,atom(M))) &
                 - vdwcut(iaL,atom(M)) - vdwsf(iaL,atom(M))*(r1 - rcut)
            Ulj = Ulj + dele
         end if
         Uel = Uel + qL*charge(M)/r1
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE

   PURE SUBROUTINE EnergyLJEL(atomfirst,atomlast,sysfirst,syslast,Uljel)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE atom_types, only: atom
      USE charges_mod, only: charge
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: Uljel
      real(wp):: rr(3),r2,r1,rr6,dele
      integer:: L,M
      Uljel = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         if (nearest_neighbor2(M,L)) cycle atom_loop
         rr(1:3) = rxyz(M,1:3) - rxyz(L,1:3)
         rr(1:2) = rr(1:2) - boxl*anint(rr(1:2)*boxli)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         r1 = sqrt(r2)
         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(atom(L),atom(M))*rr6 - bij(atom(L),atom(M))) &
                 - vdwcut(atom(L),atom(M)) - vdwsf(atom(L),atom(M))*(r1 - rcut)
            Uljel = Uljel + dele
         end if
         Uljel = Uljel + charge(L)*charge(M)/r1
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE

END MODULE lj_el_mod

!!>include 'Henrys_law_calc_mod.f90'

MODULE Henrys_law_calc_mod
    USE precision_mod, only: wp
    USE global_vars
    implicit none
    PUBLIC:: Henrys_law_calc
CONTAINS

   SUBROUTINE Henrys_law_calc(ntrial,zlower,zupper,pr,fover,Uljel,Khenry,facc)
      USE coordinates_mod
      USE atom_types
!     USE constants_mod, only: sigma_2
      USE seaton_mod, only: sigi2
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE global_vars, only: etot
      USE lj_el_mod
      USE probe_mol_mod
      integer,intent(in):: ntrial
      real(wp),intent(in):: zlower,zupper
      type(probe_mol),intent(in):: pr
      real(wp),intent(in):: fover
      real(wp),intent(out):: Uljel,Khenry,facc
      real(wp),allocatable:: probe(:,:)
      integer,allocatable:: patom(:)
      real(wp):: aa(3,3),q(4),x,y,z,rnaccept,dU,dKHenry
      integer:: it,nat
      logical:: rotate
      nat = pr%n
      allocate( probe(3,nat), patom(nat) )
      patom = pr%atom
      rotate = ((size(probe,2) > 1) .and. nat > 1)
      rnaccept = 0
      Uljel = 0.0_wp
      KHenry = 0.0_wp
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
         probe(1:2,:) = probe(1:2,:) - anint(probe(1:2,:)*boxli)*boxl
! first check if there is a large overlap
         if (overlap(fover)) cycle insert_loop
         rnaccept = rnaccept + 1.0_wp
!         call ENERGY_LJ_EL(probe,patom,nat,1-nseed,natom,dU)
         Uljel = Uljel + dU
         dKHenry = exp(-dU/etot)
         KHenry = KHenry + dKHenry
      end do insert_loop
      Uljel = Uljel/real(ntrial,wp)
      KHenry = KHenry/real(ntrial,wp)
      facc = rnaccept/real(ntrial,wp)
      RETURN
!
CONTAINS
!
      pure FUNCTION overlap(fr)
         real(wp),intent(in):: fr
         logical:: overlap
         integer:: k,jj,j,ic,nc,ncell(125),neigh
         real(wp):: dr(3)
         overlap = .false.
         probe_atom_loop: do k = 1,nat
         ic = CELL( probe(1:3,k) )
         call NEIGCELL(ic,1,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0) cycle cell_loop
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               if (atom(j) /= iSilicon) then
               dr(1:3) = probe(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               if (dot_product(dr,dr) < ( fr*( pr%rad(k) + sigi2(atom(j)) ) )**2) then
                  overlap = .true.
                  RETURN
               end if
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
         end do probe_atom_loop
      END FUNCTION overlap
   END SUBROUTINE Henrys_law_calc

END MODULE Henrys_law_calc_mod



PROGRAM porous_silica
      USE precision_mod, only: wp
      USE global_vars
      USE constants_mod
      USE atom_types
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
USE ran_point_sphere
USE quat2mat_mod
      USE hydrogens_mod
      USE frames_mod
      USE Relax_mod
      USE nlist_mod
      USE verlet_list
      USE HKNonLattice_mod
      USE deposition_mod
      USE rdf_mod
      USE Keating_parameters, only: ASiO
      USE add_water_mod
      USE relax_dangling_mod
      USE files_mod
      USE Relax_mod
      USE attachment_mod
      USE insert_seeds_mod
      USE probe_mol_mod
      USE voidage_mod
      USE surface_atoms_mod
      USE keating,only:energy4
      USE repul_energy_mod
      USE seaton_mod
      USE charges_mod
      USE lj_el_mod
      USE Henrys_law_calc_mod
      implicit none
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: aSiOSi = pi*(140.0_wp/180.0_wp)
      real(wp),parameter:: r_ox = 0.15_wp
      real(wp),parameter:: d_ox = 2.0_wp*r_ox
      real(wp):: top_atom,void_crit,rverlet,dlayer !!,z_pos
      real(wp):: dUlj,dUel,rr(3),q(4),aa(3,3),ulj(100),uel(100)
type(probe_mol):: pr1,pr2
real(wp),allocatable:: probe(:,:),PCHARGE(:)
integer,allocatable:: PATOM(:)
integer:: navg = 100
      integer:: i,L,iz,ipconfig,ntotal,imve,N_OH_groups,j,k
      real:: timep(0:10),t0,t1
      character:: date*8,time*10
      logical:: success,ok
integer:: ntrial_henry
real(wp):: uljel,khenry,facc,fover
!
      OPEN(unit=10,FILE = 'psil1.dat',STATUS = 'OLD')
!-----------------------------------------------------------
!     Read in simulation data
!     nattached = number of molecules already attached
!     ntotal = total number of molecules to be attached
!     nseed = number of initial surface oxygens
!     boxl = box side length of cell
!     ipconfig = frequency of xyz,RDF dump
!     irand_seed = random number seed
!-----------------------------------------------------------
      READ(10,*) nattached
      READ(10,*) ntotal
      READ(10,*) nseed
      READ(10,*) Etot
      READ(10,*) fwater
      READ(10,*) boxl
      READ(10,*) void_crit
      READ(10,*) del_rxn
      READ(10,*) e_activ
      READ(10,*) ipconfig
      READ(10,*) irand_seed
      READ(10,*) NRELAX
      READ(10,*) rverlet
      close(UNIT = 10,STATUS = 'KEEP')
      natom_max = 2*nseed + 7*ntotal
      allocate(rxyz(1 - nseed:natom_max,3))
      allocate(atom(1 - nseed:natom_max))
      allocate(charge(1 - nseed:natom_max))
      allocate( proximity(natom_max,4) )
      call Init_Verlet_List(rverlet,0.26_wp)
      call Init_HKNonLattice(natom_max)
      call LJ_init
!
      proximity = 0
boxl = boxl*10.0
      boxli = 1.0_wp/boxl
      boxl2 = boxl/2.0_wp
!
      OPEN(unit=15,file='psil_config.xyz')
      OPEN(unit=20,file='psil_topology.out')
      OPEN(unit=14,file='psil_RDF.out',POSITION = 'APPEND')
      OPEN(unit=21,file='psil_rand.out')
      OPEN(unit=22,file='psil_dens.out')
      OPEN(unit=23,file='psil_void_SiO4.out')
      OPEN(unit=26,file='psil_void_H2O.out')
      OPEN(unit=24,file='psil_topatom.out')
      OPEN(unit=25,file='psil_time.out')
      OPEN(unit=77,file='psil_SiO4_attempt.out')
      OPEN(unit=78,file='psil_SiO4_add.out')
      OPEN(unit=79,file='psil_H2O_add.out')
      OPEN(unit=80,file='psil_H2O_remove.out')
      OPEN(unit=81,file='psil_bond_switch.out')

!      if (nattached == 0) then
         !call insert_seeds
         !n_cluster = nseed
         !n_cluster_old = n_cluster
         !forall(i = 1:n_cluster) atomL(i) = i  ! Label each seed as a cluster
         !if (make_movie) then ! write out first frame of movie
            !call new_frame_file(imve,'frame',0)
            !call write_frame(imve,1 - nseed,natom)
         !end if
         !top_atom = 0.0_wp
      !else
         !call read_xyz(15,nattached)
         !call read_proximity(20,nattached)
         !call read_rand_state(21)
         !call HKNonLattice(natom,proximity,n_cluster,atomL)
         !n_cluster_old = n_cluster
         !write(*,*) 'n_cluster = ',n_cluster
         !write(*,*) 'atomL = '; write(*,'(10i6)') atomL(1:natom)
         !top_atom = maxval(rxyz(1:natom,3))
         !forall(i = 1:nseed) rxyz(i - nseed,:) = rxyz(i,:) - (/0.0_wp,0.0_wp,bondl_SiO4/)
         !atom(1 - nseed:0) = iSilicon
      !end if
!
      call init_probe_mols()
      call WRITE_PROBE_MOL(999,probe_SiOH4)

!      call INIT_NLIST(boxl,boxl,500.0_wp/angstrom,max(5.0_wp/angstrom,2.0_wp*maxval(probe_H2O%rad)))
      !call NEW_NLIST
      call cpu_time(timep(0))
!
!----------------------------------------------------------------------
      main_loop: do while (nattached < ntotal)
      call cpu_time(timep(1))
      nattached = nattached + 1
!
      pr1 = probe_SiOH4
      pr2 = probe_SiOH4
      natom = pr2%n
      ALLOCATE(probe(3,natom),patom(natom),pcharge(natom))

      patom = pr1%atom
      pcharge = 0.0_wp  !pr1%q
      atom(1:natom) = pr2%atom
      charge(1:natom) = 0.0_WP!pr2%q
      write(*,*) charge(1:natom)
      write(*,*) 'net charge = ',sum(charge(1:natom))

      dlayer = 1.4_wp/100
      ulj = 0.0_wp
      uel = 0.0_wp

      do k = 1,navg
      call ran4sph(q)
      call quat2mat(q,aa)
      probe = matmul(aa,probe_SiOH4%r)
      do j = 1,navg
         call ran4sph(q)
         call quat2mat(q,aa)
         pr2%r = matmul(aa,probe_SiOH4%r)
         do i = 39,100
            forall(i = 1:natom) rxyz(i,:) = pr2%r(:,i)
            rxyz(1:natom,1) = rxyz(1:natom,1) + i*dlayer
            call ENERGY_LJ_EL(probe,patom,pcharge,natom,1,natom,dUlj,dUel)
            ulj(i) = ulj(i) + dUlj
            uel(i) = uel(i) + dUel
         end do
      end do
      end do
      do i = 39,100
         write(997,*) i*dlayer,Ulj(i)/(navg*navg),Uel(i)/(navg*navg)
      end do

call cpu_time(timep(2))
PRINT *,'TIME = ',timep(2) - timep(1)
STOP




      dlayer = (5.0_wp/angstrom)
      if (use_voidage_calc) then
         call voidage_profile(probe_SiO4,dlayer)
      else
         call den_profile(1,natom,dlayer)
         voidage(1:nbin) = voidage_approx(nden(1:nbin,0),2.0_wp*r_tet*angstrom)
      end if
      call prob_dist(void_crit)
      itop_nonp_SiO4 = itop_nonp
!
! call Henrys_law_calc(ntrial,zlower,zupper,pr,fover,Uljel,Khenry,facc)
call cpu_time(t0)
ntrial_henry = 10000
fover = 0.8_wp
 call Henrys_law_calc(ntrial_henry,0.5_wp,1.0_wp,probe_SiOH4,fover,Uljel,Khenry,facc)
call cpu_time(t1)
print *
print *, ' fover = ',fover
print *, ' Uljel = ',Uljel
print *, 'Khenry = ',Khenry
print *, '  facc = ',facc
print *, ntrial_henry,fover,natom,t1 - t0

call cpu_time(t0)
ntrial_henry = 10000
fover = 0.95_wp
 call Henrys_law_calc(ntrial_henry,0.5_wp,1.0_wp,probe_SiOH4,fover,Uljel,Khenry,facc)
call cpu_time(t1)
print *
print *, ' fover = ',fover
print *, ' Uljel = ',Uljel
print *, 'Khenry = ',Khenry
print *, '  facc = ',facc
print *, ntrial_henry,fover,natom,t1 - t0
call cpu_time(t0)
ntrial_henry = 1000
fover = 0.8_wp
 call Henrys_law_calc(ntrial_henry,1.0_wp,1.5_wp,probe_SiOH4,fover,Uljel,Khenry,facc)
call cpu_time(t1)
print *
print *, ' Uljel = ',Uljel
print *, 'Khenry = ',Khenry
print *, '  facc = ',facc
print *, ntrial_henry,fover,natom,t1 - t0


do i = 1,int(maxval(rxyz(1:natom,3))/dlayer) + 1
call Henrys_law_calc(ntrial_henry,(i - 1)*dlayer,i*delz,probe_SiOH4,fover,Uljel,Khenry,facc)
write(*,*) (i - 0.5_wp)*delz,Uljel,Khenry,facc
end do




!
      do ! attach SiO3
         do
            call sample_z_bin(iz)
            !! call sample_z_bin_uniform(iz)
            call get_atom_for_attachment2(iz,d_ox,aSiO,L)
            if (L > 0) exit
            print *,'L == 0',iz
         end do
         write(77,*) nattached,iz,rxyz(L,3)
         call attach_tetrahedra(L,aSiOSi,aSiO,success)
         if (success) exit
         !!write(*,*) 'failed to attach SiO3'
      end do
call adjust_hydrogens()
      write(78,*) nattached,iz,rxyz(L,3)
      call cpu_time(timep(2))

      call new_vlist
!
!     call RELAX(NRELAX, nseed+1, natom, hot_atoms=.TRUE.)
!!print *,'energy before RELAX ' !,energy4(1,natom) + repul_energy(1,natom)
      call RELAX(NRELAX, nseed + 1, natom, delt = 0.5_wp, &
                 alist = (/natom - 3,natom - 2,natom - 1,natom/),na = 4, &
                 quench = .FALSE., scale_repul = .FALSE.)
!!print *,'energy after RELAX  ' !,energy4(1,natom) + repul_energy(1,natom)
!
!     Local relax
      call relax_dangling(3)

      call cpu_time(timep(3))

!!print *,'energy after LOCAL RELAX  ' !,energy4(1,natom) + repul_energy(1,natom)
      N_OH_groups = count(atom(1:natom) == iOxygenH)
!     Global relax
      call relax_dangling(N_OH_groups)
call adjust_hydrogens()

!!print *,'energy after GLOBAL RELAX  ' !,energy4(1,natom) + repul_energy(1,natom)
      call cpu_time(timep(4))
! try to add water
      call GET_NUM_ADD_WATER_ATTEMPTS()
      do i = 1,ntry_add_water
!!print *,'trying to add water',i
call assign_charge(1,natom)

         if (use_voidage_calc) then
            call voidage_profile(probe_H2O,dlayer)
         else
            call den_profile(1,natom,dlayer)
            voidage(1:nbin) = voidage_approx(nden(1:nbin,0),2.0_wp*r_ox*angstrom)
         end if
         call prob_dist(void_crit)
         call sample_z_bin(iz)
         call add_water(iz,success)
         if (success) then
            N_OH_groups = count(atom(1:natom) == iOxygenH)
            ! Global relax
            call relax_dangling(N_OH_groups)
!!print *,'energy after  add_water'!,energy4(1,natom) + repul_energy(1,natom)
         end if
      end do

call adjust_hydrogens()
      call cpu_time(timep(5))
      call date_and_time(date,time)

!     write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
!     write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
      write(* ,'(i7,a10,2x,a10,2(3x,f0.3))') natom,date,time,timep(2) - timep(1), &
         timep(3) - timep(2),timep(4) - timep(3),timep(5) - timep(4),timep(5) - timep(0)
      write(25,'(i7,a10,2x,a10,6(2x,f0.3))') natom,date,time,timep(2) - timep(1), &
         timep(3) - timep(2),timep(4) - timep(3),timep(5) - timep(4),timep(5) - timep(0)
      top_atom = maxval(rxyz(1:natom,3))
      write(24,'(2i7,3f12.6)') nattached,natom,top_atom
!
      if (make_movie) then
         call new_frame_file(imve,'frame',nattached)
         call write_frame(imve,1,natom)
      end if
!
      if ( mod(nattached,ipconfig) == 0 .or. (nattached == ntotal) ) then
!--------Periodically check the connectivity and write out xyz coordinates,
!--------RDF, Connectivity, etc.
         call check_proximity(ok,i)
         if (.NOT.ok) then
            write(*,*) 'proximity array is NOT consistent'
            write(*,*) 'check_proximity ',ok,i
            stop
         end if
         call check_proximity2(ok,i)
         if (.NOT.ok) then
           write(*,*) 'proximity array is NOT consistent'
           write(*,*) 'check_proximity2 ',ok,i
           stop
         end if
         if (ANY(atom(proximity(1:nseed,1)) == iOxygen)) then
            write(*,*) 'seeds should not be bonded to Oxygens'
            stop
         end if
         call rdf_calc
         call rdf_print(14,nattached)
         call write_xyz(15,nattached)
         rewind(15)
         call write_proximity(20,nattached)
         rewind(20)
         call write_rand_state(21)
         rewind(21)
         if (use_voidage_calc) call den_profile(1,natom,(5.0_wp/angstrom))
         call write_den_profile(22)
         call voidage_profile(probe_SiO4,dlayer)
         call write_voidage(23)
         call SURFACE_ATOMS(nseed + 1,natom,(1.5_wp/angstrom))
         write(*,*) 'fraction of surface atoms = ',fsurfa
         !call write_surface_atoms_zbin(191)
         call voidage_profile(probe_H2O,dlayer)
         call write_voidage(26)
      end if

      end do main_loop
!----------------------------------------------------------------------
      close(14)
      close(15)
      close(20)
      close(21)
      close(22)
      close(23)
      STOP
END PROGRAM Porous_silica

