
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

!!>include 'global_vars.f90'

MODULE global_vars
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, kboltzT
   integer:: natom,nseed,natom_max
   integer:: nattached, nrelax
END MODULE

!!>include 'rand_mod.f90'

MODULE rand_mod
! a random number generator (Numerical Recipies)
! with auxillary subroutines for storing and setting the
! state space
   use precision_mod, only: wp,i4b
   implicit none
   private
   public:: rand,get_rand_state,set_rand_state,read_rand_state,write_rand_state
   integer(i4b),private,save:: ix = -1, iy = -1
   integer(i4b),public,save:: irand_seed = -1
!
   contains

   subroutine set_rand_state(ix0,iy0)
      integer(i4b),intent(in):: ix0,iy0
      real(wp):: tmp
      tmp = rand()
      ix = ix0
      iy = iy0
   end subroutine set_rand_state

   subroutine get_rand_state(ix0,iy0)
      integer(i4b),intent(out):: ix0,iy0
      ix0 = ix
      iy0 = iy
   end subroutine get_rand_state

   subroutine read_rand_state(iu)
      integer(i4b),intent(in):: iu
      real(wp):: tmp
      tmp = rand()
      read(iu,*) ix,iy
      rewind(iu)
   end subroutine read_rand_state

   subroutine write_rand_state(iu)
      integer(i4b),intent(in):: iu
      write(iu,*) ix,iy
      call flush(iu)
   end subroutine write_rand_state

   function rand()
!-----random number generator
      real(wp):: rand
      integer(i4b),parameter::ia=16807,im=2147483647,iq=127773,ir=2836
      integer(i4b),save:: k
      real(wp),save::am
      if(irand_seed <= 0 .or. iy < 0 )then
         am = nearest(1.0_wp,-1.0_wp)/im
         iy = ior(ieor(888889999,abs(irand_seed)),1)
         ix = ieor(777755555,abs(irand_seed))
         irand_seed = abs(irand_seed)+1
      end if
      ix = ieor(ix,ishft(ix,13))
      ix = ieor(ix,ishft(ix,-17))
      ix = ieor(ix,ishft(ix,5))
      k = iy/iq
      iy = ia*(iy-k*iq)-ir*k
      if(iy < 0)iy = iy+im
      rand = am*ior(iand(im,ieor(ix,iy)),1)
   end function rand

! 'minimal standard' RNG
!   function rand0()
!      real(wp):: rand0
!      real(wp),parameter:: im=2147483647.0_wp,ia=16807.0_wp
!      real(wp),parameter:: am=1.0_wp/im
!      irand_seed = mod(ia*mod(ia*irand_seed,im),im)
!      rand0 = am*irand_seed
!   end function

END MODULE rand_mod

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
1     if(ir-L < M)then
        do j=L+1,ir
          a=arr(j)
          do i=j-1,1,-1
            if(arr(i) <= a)goto 2
            arr(i+1)=arr(i)
          end do
          i=0
2         arr(i+1)=a
        end do
        if(jstack == 0)return
        ir=istack(jstack)
        L=istack(jstack-1)
        jstack=jstack-2
      else
        k=(L+ir)/2
        temp=arr(k)
        arr(k)=arr(L+1)
        arr(L+1)=temp
        if(arr(L+1) > arr(ir))then
          temp=arr(L+1)
          arr(L+1)=arr(ir)
          arr(ir)=temp
        end if
        if(arr(L) > arr(ir))then
          temp=arr(L)
          arr(L)=arr(ir)
          arr(ir)=temp
        end if
        if(arr(L+1) > arr(L))then
          temp=arr(L+1)
          arr(L+1)=arr(L)
          arr(L)=temp
        end if
        i=L+1
        j=ir
        a=arr(L)
3       continue
          i=i+1
        if(arr(i) < a)goto 3
4       continue
          j=j-1
        if(arr(j) > a)goto 4
        if(j < i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(L)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack > NSTACK) stop 'NSTACK too small in qsort'
        if(ir-i+1 >= j-L)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=L
          L=i
        end if
      end if
      goto 1
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
      if(inc <= n) goto 1
2     continue   ! Loop over the partial sorts
      inc = inc/3
      do i=inc+1,n   ! Outer loop of straight insertion.
         b=v(i)
         j=i
3        if(v(j-inc) > b)then   ! Inner loop of straight insertion.
            v(j)=v(j-inc)
            j=j-inc
            if(j <= inc)goto 4
            goto 3
         endif
4        v(j)=b
      enddo
      if(inc > 1)goto 2
      return
   END SUBROUTINE

   SUBROUTINE sort3(iv)
      integer,intent(inout):: iv(:)
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
      if (iv(3) < iv(2)) call swap(iv(3),iv(2))
      if (iv(2) < iv(1)) call swap(iv(2),iv(1))
   contains
      subroutine swap(x, y)
         integer,intent(inout):: x,y
         integer:: tmp
         tmp = x
         x = y
         y = tmp
      end subroutine
   END SUBROUTINE

END MODULE

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

!!>include 'seaton_mod.f90'

MODULE seaton_mod
   USE precision_mod
   USE global_vars, only: angstrom, K_ev, qstar, pi
   implicit none
! atom types
   integer,parameter:: iO_Sil=0
   integer,parameter:: iSi_Sil=1
   integer,parameter:: iO_OH_sil=2
   integer,parameter:: iH_OH_sil=3
   integer,parameter:: iO_H2O=4
   integer,parameter:: iH_H2O=5
   integer,parameter:: iO_CO2=6
   integer,parameter:: iC_CO2=7
   integer,parameter:: iN_N2=8
   integer,parameter:: iO_O2=9
   integer,parameter:: iN_N2charge=10
   integer,parameter:: iO_O2charge=11
   integer,parameter:: ntyplj = 11  ! number of LJ types
   real(wp),parameter:: eps_O_sil   =185.0_wp*K_ev, sig_O_sil  = 2.708_wp/angstrom,q_O_sil = -0.64025_wp*qstar
   real(wp),parameter:: eps_Si_sil  =  0.0_wp*K_ev, sig_Si_sil = 0.0_wp/angstrom,  q_Si_sil= -2.0_wp*q_O_sil
   real(wp),parameter:: eps_O_OH_sil=185.0_wp*K_ev, sig_O_OH_sil= 3.0_wp/angstrom, q_O_OH_sil= -0.533_wp*qstar
   real(wp),parameter:: eps_H_OH_sil=  0.0_wp*K_ev, sig_H_OH_sil= 0.0_wp/angstrom, q_H_OH_sil=  0.206_wp*qstar
   real(wp),parameter:: eps_H_H2O   =  0.0_wp*K_ev,  sig_H_H2O = 0.0_wp/angstrom,   q_H_H2O=  0.417_wp*qstar
   real(wp),parameter:: eps_O_H2O   = 76.58_wp*K_ev, sig_O_H2O = 3.1506_wp/angstrom,q_O_H2O= -2.0_wp*q_H_H2O
   real(wp),parameter:: eps_O_CO2   = 82.997_wp*K_ev,sig_O_CO2 = 3.064_wp/angstrom, q_O_CO2= -0.33225_wp*qstar
   real(wp),parameter:: eps_C_CO2   = 29.999_wp*K_ev,sig_C_CO2 = 2.785_wp/angstrom, q_C_CO2= -2.0_wp*q_O_CO2
   real(wp),parameter:: eps_N_N2    = 34.897_wp*K_ev, sig_N_N2 = 3.3211_wp/angstrom, q_N_N2= -0.5475_wp*qstar
   real(wp),parameter:: eps_O_O2    = 43.183_wp*K_ev, sig_O_O2 = 3.1062_wp/angstrom, q_O_O2= -0.3577_wp*qstar
!  bond lengths
   real(wp),parameter:: bondl_SiO = 1.60_wp/angstrom
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
! LJ parameters
   real(wp),parameter:: epsi(0:ntyplj)=(/ &
                        eps_O_sil, &
                        eps_Si_sil, &
                        eps_O_OH_sil, &
                        eps_H_OH_sil, &
                        eps_O_H2O, &
                        eps_H_H2O, &
                        eps_O_CO2, &
                        eps_C_CO2, &
                        eps_N_N2, &
                        eps_O_O2, &
                        0.0_wp, &
                        0.0_wp /)
   real(wp),parameter:: sigi(0:ntyplj)=(/ &
                        sig_O_sil, &
                        sig_Si_sil, &
                        sig_O_OH_sil, &
                        sig_H_OH_sil, &
                        sig_O_H2O, &
                        sig_H_H2O, &
                        sig_O_CO2, &
                        sig_C_CO2, &
                        sig_N_N2, &
                        sig_O_O2, &
                        0.0_wp, &
                        0.0_wp /)
   real(wp),parameter:: sigi2(0:ntyplj)=sigi*0.5_wp
! Charges
   real(wp),parameter:: q_seaton(0:11)=(/ &
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

!!>include 'constants_mod.f90'

MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars, only: angstrom,pi
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
   integer:: n_O2,n_N2,n_CO2
   real(wp):: energy
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
   real(wp),parameter:: sig_He   = 2.64_wp/angstrom
! the following are now taken from Seaton
!   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
!   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
!   real(wp),parameter:: sigma_Obridge = 2.7_wp/angstrom
!   real(wp),parameter:: sigma_OH = 3.00_wp/angstrom
!   real(wp),parameter:: sigma_O = 2.94_wp/angstrom
!   real(wp),parameter:: sigma_Si = 4.20_wp/angstrom
!   real(wp),parameter:: sigma_C = 2.70_wp/angstrom
!   real(wp),parameter:: sigma_N = 3.296_wp/angstrom
!   real(wp),parameter:: bondl_O2 = 1.169_wp/angstrom
!   real(wp),parameter:: bondl_N2 = 1.097_wp/angstrom
!   real(wp),parameter:: bondl_SiO = 1.62_wp/angstrom
!   real(wp),parameter:: bondl_CO2 = 0.958_wp/angstrom
! united atom SiOH4
   real(wp),parameter:: r_tet = bondl_SiO + sigma_OH*0.5_wp
!
   real(wp),parameter:: sigma(0:3) = (/ sig_O_sil, sig_Si_sil,sigma_OH, sig_H_OH_sil /)
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
   integer,parameter:: iOw=4
   integer,parameter:: iHw=5
   integer,parameter:: ntyp = 5
   character(2),parameter:: atom_name(0:ntyp)=(/ 'O ','Si','OH','H ','Ow','Hw' /)
   integer,parameter:: ncmax(0:ntyp)=(/2,4,2,1,1,2/)
   integer,allocatable:: atom(:),Ox_atom(:,:),N2_atom(:,:),CO2_atom(:,:)
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

END MODULE

!!>include 'coordinates_mod.f90'


MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:),Ox_xyz(:,:,:)
   real(wp),allocatable:: Ox_fxyz(:,:,:),N2_fxyz(:,:,:),CO2_fxyz(:,:,:)
   real(wp),allocatable:: Imige_Ox_xyz(:,:,:),N2_xyz(:,:,:),CO2_xyz(:,:,:)
   real(wp),allocatable:: Imige_CO2_xyz(:,:,:),Imige_N2_xyz(:,:,:)
   real(wp),allocatable:: Ox_vxyz(:,:,:),N2_vxyz(:,:,:),CO2_vxyz(:,:,:),mol_vec(:,:)

   real(wp):: boxl,boxli,boxl2
   interface pbc
      module procedure  pbc_v, pbc_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      r(1:3) = r(1:3) - boxl*anint(r(1:3)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      r(:,1:3) = r(:,1:3) - boxl*anint(r(:,1:3)*boxli)
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
      if (it /= nattached)stop 'it /= nattached'
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
               return
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
               if (atom(k)==iSilicon) GOTO 100   ! but not Si
               if (atom(k)==iHydrogen) GOTO 100  ! and not H
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

!!>include 'atom_list_mod.f90'

MODULE atom_list_mod
!
   integer,private,parameter:: n_atoml_max = 4000
   TYPE atom_list
      integer:: i(n_atoml_max)
      integer:: n = 0
   END TYPE
!
CONTAINS

   SUBROUTINE add_to_list(lst,iat)
      type(atom_list),intent(inout):: lst
      integer,intent(in):: iat
      lst%n = lst%n + 1
      lst%i(lst%n) = iat
   END SUBROUTINE

   SUBROUTINE remove_from_list(lst,iat)
      type(atom_list),intent(inout):: lst
      integer,intent(in):: iat
      integer:: i
      do i = 1,lst%n
         if(lst%i(i)==iat)exit
      end do
      lst%i(i) = lst%i(lst%n)
      lst%i(lst%n) = 0
      lst%n = lst%n - 1
   END SUBROUTINE

   SUBROUTINE rand_from_list(lst,iat)
      USE rand_mod
      integer,intent(out):: iat
      type(atom_list),intent(in):: lst
      integer:: j
      j = int(rand()*lst%n) + 1
      iat = lst%i(j)
   END SUBROUTINE

END MODULE atom_list_mod

!!>include 'charges2_mod.f90'

MODULE charges_mod
   USE precision_mod
   USE global_vars, only: nseed
   USE atom_types
   USE seaton_mod, only: qi => q_seaton
   implicit none
   real(wp),allocatable:: charge(:),Ox_gas_charge(:,:),N2_gas_charge(:,:),CO2_gas_charge(:,:)
   PUBLIC:: charge, assign_charge, qsi
!
CONTAINS

   pure function qsi(n)
      real(wp):: qsi
      integer,intent(in):: n
      qsi = (n-4)*qi(iOxygen)*0.5_wp - n*(qi(iOxygenH) + qi(iHydrogen))
   end function

   SUBROUTINE assign_charge(ifirst,ilast)
      USE connectivity_mod, only: noh
      USE atom_types
      integer,intent(in):: ifirst,ilast
      real(wp):: qo,qh,qoh
      integer:: i
      qoh= qi(iOxygenH) + qi(iHydrogen)
      qh = qi(iHydrogen)
      qo = qi(iOxygen)
      do i = ifirst,ilast
         select case(atom(i))
         case(iSilicon)
            charge(i) = qsi(noh(i))
         case(iOxygen)
            charge(i) = qo
         case(iOxygenH)
            charge(i) = qoh
         case(iHydrogen)
            STOP 'no explicit H used'
            !charge(i) = qh
         case default
            stop 'assign_charge: unknown atom type'
         end select
      end do
   END SUBROUTINE

END MODULE charges_mod

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
         a12 = -a*b/(1.0_wp + c)    ! (a*b*(c-1))/(a*a + b*b)
         aa(1,1:3) = (/ a11, a12, a /)
         aa(2,1:3) = (/ a12, a22, b /)
         aa(3,1:3) = (/  -a,  -b, c /)
      else  ! r1  is close to +/- (0,0,1)
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

!!>include 'nlist_mod.f90'

MODULE NLIST_MOD
      USE precision_mod, only: wp
      USE global_vars, only: natom_max,natom
      USE coordinates_mod, only: rxyz,boxl2
      implicit none
      integer,public:: neigh
      real(wp):: delx,dely,delz,delmin
      real(wp),private:: delxi,delyi,delzi
      integer:: ncelx,ncely,ncelz,ncelt
      integer,parameter:: neighmx=350
      integer,allocatable,public:: hoc(:),hoc_old(:)
      integer,public:: ncell(neighmx)
      integer,allocatable,public:: ll_old(:),ll(:),lr(:)
      public:: NEIGCELL,cell,INIT_NLIST,NEW_NLIST,print_cell
      public:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll
!
!     ll(i)      : linked list particle i
!     hoc(ic)    : head of chain cell ic
!     ncelx,y,z  : number of cells in x, y or z direction
!     ncelt      : total number of cells
!     neigh      : number of cells for interactions
!
CONTAINS

   pure function nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   end function

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
      write(*,*)'ncelx = ',ncelx
      write(*,*)'ncely = ',ncely
      write(*,*)'ncelz = ',ncelz
      write(*,*)'ncelt = ',ncelt
      write(*,*)'delx = ',delx
      write(*,*)'dely = ',dely
      write(*,*)'delz = ',delz
      write(*,*)'delmin = ',delmin
      allocate(hoc(ncelt),hoc_old(ncelt))
      LL = 0
      HOC = 0
      LR = 0
   END SUBROUTINE

   PURE FUNCTION CELL(XYZ)
!     determines cell number for position x,y,z
      real(wp),intent(in)::  XYZ(:)
      integer:: CELL
      CELL = INT((XYZ(1)+boxl2)*delxi) &
           + INT((XYZ(2)+boxl2)*delyi)*ncelx &
           + INT((XYZ(3)+boxl2)*delzi)*ncelx*ncely + 1
      RETURN
   END FUNCTION CELL

   SUBROUTINE NEW_NLIST
!     makes a new neighbour list using the linked-list algorithm
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head-of-chain
      ! make linked list:
      do i = 1,natom
         ! determine cell number
         ic = CELL(rxyz(i,:))
if(ic < 1) stop 'NEW_NLIST: ic < 1'
if(ic > ncelt) stop 'NEW_NLIST: ic > ncelt'
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
            iccz = iccz + ncelz
         else if (iccz > ncelz) then
            iccz = iccz - ncelz
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
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp,-bondl /)
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
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp,-bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

   SUBROUTINE SiOH4_coords(bondl,bondl_OH,angle_SiOH,tetra)
      USE rotate_axis_mod
      USE rand_mod
      real(wp),intent(in):: bondl,bondl_OH,angle_SiOH
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
      cang = cos(angle_SiOH)
      sang = sin(angle_SiOH)
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp,-bondl /)
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
         tetra(1:3,j+4) = tetra(1:3,j) + matmul(aa,(/ xx,yy,zz /))
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
      TQ = TY+TZ
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
         qx = (2.0_wp*rand()-1.0_wp)
         qy = (2.0_wp*rand()-1.0_wp)
         qz = (2.0_wp*rand()-1.0_wp)
         s1=qx*qx+qy*qy+qz*qz
         if (s1 <= 1.0_wp)exit
      end do
      return
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
      return
   END SUBROUTINE

   SUBROUTINE ran3sph2(x,y,z)
      real(wp),parameter::twopi=6.283185307179586476925286766559005768394_wp
      real(wp),intent(out)::x,y,z
      real(wp)::s,tmp
      z = 2.0_wp*rand()-1.0_wp
      tmp = rand()*twopi
      s = sqrt(1.0_wp-z*z)
      x = s*cos(tmp)
      y = s*sin(tmp)
      return
   END SUBROUTINE

   SUBROUTINE ran_point_on_sphere(qx,qy,qz)
      real(wp),intent(out):: qx,qy,qz
      real(wp)::s1,r1
      do
        qx = (2.0_wp*rand()-1.0_wp)
        qy = (2.0_wp*rand()-1.0_wp)
        qz = (2.0_wp*rand()-1.0_wp)
        s1=qx*qx+qy*qy+qz*qz
        if (s1 <= 1.0_wp)exit
      end do
      r1 = 1.0_wp/sqrt(s1)
      qx = qx*r1
      qx = qx*r1
      qz = qx*r1
      return
   END SUBROUTINE

   SUBROUTINE ran4sph(q)
! random point on surface of 4-Sphere (a Quaternion)
      real(wp),intent(out)::q(4)
      real(wp)::s1,s2
      do
        q(1) = (2.0_wp*rand()-1.0_wp)
        q(2) = (2.0_wp*rand()-1.0_wp)
        s1 = q(1)*q(1)+q(2)*q(2)
        if (s1 <= 1.0_wp)exit
      end do
      do
        q(3) = (2.0_wp*rand()-1.0_wp)
        q(4) = (2.0_wp*rand()-1.0_wp)
        s2 = q(3)*q(3) + q(4)*q(4)
        if (s2 <= 1.0_wp)exit
      end do
      q(3) = q(3)*sqrt((1.0_wp-s1)/s2)
      q(4) = q(4)*sqrt((1.0_wp-s1)/s2)
      return
   END SUBROUTINE

END MODULE

!!>include 'probe_mol_mod_DVF.f90'

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
      real(wp):: r(3,9)
      integer:: atom(9)
      real(wp):: rad(9),q(9)
   END TYPE probe_mol
!
   type(probe_mol):: probe_SiOH4, probe_SiO4, probe_tip3p,probe_He
   type(probe_mol):: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet
!
CONTAINS

   SUBROUTINE new_probe_mol(P,n,radius,atomtype)
      type(probe_mol),intent(out):: P
      integer,intent(in):: n
      real(wp),intent(in),optional:: radius
      integer,intent(in),optional:: atomtype
      P%n = n
!      allocate(P%r(3,n),P%atom(n),P%rad(n),P%q(n))
      P%r = 0.0_wp
      P%rad = 0.0_wp
      P%atom = 0
      if(present(radius)) P%rad = radius
      if(present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types
      USE seaton_mod
      USE charges_mod
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigma_H2O*0.5_wp)
! Helium
      call new_probe_mol(probe_He,1,sig_He*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom(1:3) = (/ iO_H2O, iH_H2O, iH_H2O /)
      probe_TIP3P%rad(1:3) = sigi2(probe_TIP3P%atom(1:3))
      probe_TIP3P%q(1:3) = q_seaton(probe_TIP3P%atom(1:3))
! Si(OH)4
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO,bondl_OH,angle_SiOH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigi2(probe_SiOH4%atom)
      probe_SiOH4%q = q_seaton(probe_SiOH4%atom)
      probe_SiOH4%q(1) = qsi(4)
! SiO4, just the 4 OH Oxygens
      call new_probe_mol(probe_SiO4,4,sigma_OH*0.5_wp)
      call tetra_coords(bondl_SiO,probe_SiO4%r)
      probe_SiO4%atom(1:4) = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad(1:4) = sigi2(probe_SiO4%atom(1:4))
! O2
      call new_probe_mol(probe_O2,3,sig_O_O2*0.5_wp)
      probe_O2%atom(1:3) = (/ iO_O2, iO_O2, iO_O2charge /)
      probe_O2%rad(3) = 0.0_wp
      probe_O2%q(1:3) = q_seaton(probe_O2%atom(1:3))
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_O2*0.5_wp /)
      probe_O2%r(:,3) = (/ 0.0_wp,0.0_wp,0.0_wp /) 
! N2
      call new_probe_mol(probe_N2,3,sig_N_N2*0.5_wp)
      probe_N2%atom(1:3) = (/ iN_N2, iN_N2, iN_N2charge /)
      probe_N2%rad(3) = 0.0_wp
      probe_N2%q(1:3) = q_seaton(probe_N2%atom(1:3))
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_N2*0.5_wp /)
      probe_N2%r(:,3) = (/ 0.0_wp,0.0_wp,0.0_wp /)
! CO2
      call new_probe_mol(probe_CO2,3,sig_O_CO2*0.5_wp)
      probe_CO2%rad(1) = sig_C_CO2*0.5_wp
      probe_CO2%atom(1:3) = (/ iC_CO2, iO_CO2, iO_CO2 /)
      probe_CO2%q(1:3) = (/ q_C_CO2,q_O_CO2,q_O_CO2 /)
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp,-bondl_CO2 /)
!
!! Si(OH)4
!      write(*,*) 'SiOH4'
!      write(*,*) probe_SiOH4%n
!      write(*,'(9i6)') probe_SiOH4%atom
!      write(*,'(9f12.6)') probe_SiOH4%rad
!      write(*,'(9f12.6)') probe_SiOH4%q
!      write(*,'(3f12.6)') probe_SiOH4%r
!! H2O
!      write(*,*) 'H2O TIP3P'
!      write(*,*) probe_TIP3P%n
!      write(*,'(9i6)') probe_TIP3P%atom
!      write(*,'(9f12.6)') probe_TIP3P%rad
!      write(*,'(9f12.6)') probe_TIP3P%q
!      write(*,'(3f12.6)') probe_TIP3P%r
!! SiO4
!      write(*,*) 'SiO4'
!      write(*,*) probe_SiO4%n
!      write(*,'(9i6)') probe_SiO4%atom
!      write(*,'(9f12.6)') probe_SiO4%rad
!      write(*,'(9f12.6)') probe_SiO4%q
!      write(*,'(3f12.6)') probe_SiO4%r
! O2
      write(*,*) 'O2'
      write(*,*) probe_O2%n
      write(*,'(9i6)') probe_O2%atom(1:3)
      write(*,'(9f12.6)') probe_O2%rad(1:3)
      write(*,'(9f12.6)') probe_O2%q(1:3)
      write(*,'(3f12.6)') probe_O2%r(1:3,1:3)
! N2
      write(*,*) 'N2'
      write(*,*) probe_N2%n
      write(*,'(9i6)') probe_N2%atom(1:3)
      write(*,'(9f12.6)') probe_N2%rad(1:3)
      write(*,'(9f12.6)') probe_N2%q(1:3)
      write(*,'(3f12.6)') probe_N2%r(1:3,1:3)
! CO2
      write(*,*) 'CO2'
      write(*,*) probe_CO2%n
      write(*,'(9i6)') probe_CO2%atom(1:3)
      write(*,'(9f12.6)') probe_CO2%rad(1:3)
      write(*,'(9f12.6)') probe_CO2%q(1:3)
      write(*,'(3f12.6)') probe_CO2%r(1:3,1:3)
   END SUBROUTINE init_probe_mols

   SUBROUTINE h2o_coords(bondl,ang,r)
      real(wp),intent(in):: bondl,ang
      real(wp),intent(out):: r(3,3)
! Coordinates for water molecule
      r(1:3,1) = 0.0_wp
      r(1:3,2) = (/ 0.0_wp, 0.0_wp, bondl /)
      r(1:3,3) = (/ 0.0_wp, bondl*sin(ang), bondl*cos(ang) /)
   END SUBROUTINE

   SUBROUTINE WRITE_PROBE_MOL(iu,P)
      USE atom_types
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      integer:: i
      write(iu,*) P%n
      write(iu,*) 'Molecule'
      do i = 1,P%n
!        write(iu,'(a6,i5,a4,2x,a3,i6,4x,3f8.3')'HETATM',i,atom_name(P%atom(i)), &
!                                               '   ',i,P%r(1:3,i)*angstrom
         write(iu,'(a2,3(1x,f14.8))') atom_name(P%atom(i)),P%r(1:3,i)*Angstrom
      end do
   END SUBROUTINE

END MODULE probe_mol_mod


!!>include 'voidage4_mod.f90'

MODULE voidage_mod
   USE precision_mod
   implicit none
CONTAINS

   FUNCTION voidage_calc(ntrial,xl,xu,yl,yu,zl,zu,pr)
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
      real(wp),intent(in):: xl,xu,yl,yu,zl,zu
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
         x = xl + rand()*(xu - xl)
         y = yl + rand()*(yu - yl)
         z = zl + rand()*(zu - zl)
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
         if (int(rnaccept) < 1000) write(30,"(a,3g14.6)")'He ',x*angstrom,y*angstrom,z*angstrom
      end do insert_loop
      voidage_calc = rnaccept/real(ntrial,wp)
   END FUNCTION voidage_calc

END MODULE voidage_mod

!!>include 'deposition4_mod_H.f90'

MODULE deposition_mod
   !USE Henrys_law_calc_mod
   USE precision_mod, only: wp
   implicit none
   logical,parameter:: use_voidage_calc = .TRUE.
   integer,parameter,private:: ntrial = 10000
   integer,parameter,private:: nbinmax = 14
   integer,parameter:: iOSiOH3 = 6
   integer,parameter:: iOSiO = 7
   integer,parameter:: iOHB = 8
   integer:: nbin,itop_nonp,itop_nonp_SiO4,itop_nonp_H2O
   real(wp):: nden(nbinmax,0:4),voidage(nbinmax,nbinmax,nbinmax)
   real(wp):: ndenrxn(nbinmax,0:8)
   real(wp):: henryk(nbinmax,nbinmax,nbinmax,4),inthenryk(nbinmax,4),volbin
   real(wp):: Eljel(nbinmax,nbinmax,nbinmax,4)
   real(wp):: intvoidage(nbinmax)
CONTAINS

   elemental function voidage_approx(no,sigmap)
      real(wp):: voidage_approx
      real(wp),intent(in):: no,sigmap
      real(wp),parameter:: alpha = 34.7_wp
      real(wp),parameter:: ro = 1.5_wp
      real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
      real(wp):: t1,t2
      t1 = -(4.0_wp/3.0_wp)*pi*no*(ro*ro*ro)
      t2 = 1.0_wp + alpha*no*(sigmap/(2.0_wp*ro))
      voidage_approx = exp(t1*t2*t2*t2)
   end function voidage_approx


   SUBROUTINE den_profile(ibegin,iend,dl)
      USE coordinates_mod
      USE atom_types
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,delzi,delz
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


   SUBROUTINE integrate_fn(ibegin,iend,fv,ifv)
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
      ifv(ibegin:iend)=ifv(ibegin:iend)/tot
   END SUBROUTINE


   SUBROUTINE voidage_profile(probe,ntrial)
      USE voidage_mod
      USE coordinates_mod
      USE probe_mol_mod
      USE nlist_mod
      type(probe_mol),intent(in):: probe
      integer,intent(in):: ntrial
      real(wp):: xl,xu,yl,yu,zl,zu
      integer:: i,j,k,nbx,nby,nbz
real::t0,t1
      nbx = ncelx
      nby = ncely
      nbz = ncelz
      do k = 1,nbz
         zl = (k-1)*delz - boxl2
         zu = k*delz - boxl2
         do j = 1,nby
            yl = (j-1)*dely - boxl2
            yu = j*dely - boxl2
            do i = 1,nbx
               xl = (i-1)*delx - boxl2
               xu = i*delx - boxl2
               voidage(i,j,k) = voidage_calc(ntrial,xl,xu,yl,yu,zl,zu,probe)
            end do
         end do
print '(6i5)',k,j-1,i-1,size(voidage,1),size(voidage,2),size(voidage,3)
!print '(3i5,4g16.8)',k,j-1,i-1,voidage(i,j,k)    
      end do
   END SUBROUTINE


   SUBROUTINE prob_dist(fv,intfv,void_crit,itop_nonp)
      real(wp),intent(in):: fv(:)
      real(wp),intent(inout)::intfv(:)
      real(wp),intent(in):: void_crit
      integer,intent(out):: itop_nonp
      integer:: i
      itop_nonp = 1
      do i = nbin,1,-1
         if ( fv(i) < void_crit ) then
            itop_nonp = i+1
            exit
         end if
      end do
      call integrate_fn(itop_nonp,nbin,fv,intfv)
   END SUBROUTINE


   SUBROUTINE sample_z_bin(itop_nonp,iz)
      USE rand_mod
      integer,intent(in):: itop_nonp
      integer,intent(out):: iz
      integer:: iv1(1)
      real(wp):: rnd
      rnd = rand()
      iv1 = MINLOC(ABS(intvoidage(itop_nonp:nbin)-rnd))+itop_nonp-1
      iz = iv1(1)
   END SUBROUTINE

   SUBROUTINE sample_z_bin_uniform(itop_nonp,iz)
      USE rand_mod
      integer,intent(in):: itop_nonp
      integer,intent(out):: iz
      real(wp):: rnd
      rnd = rand()
      iz = int(rnd*(nbin - itop_nonp + 1)) + itop_nonp
   END SUBROUTINE

END MODULE deposition_mod



PROGRAM SIMBOX
      USE precision_mod
      USE global_vars
      USE constants_mod
      USE atom_types
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      USE files_mod
      USE probe_mol_mod
      USE seaton_mod
      USE rand_mod
      USE atom_types
      USE nlist_mod
      USE charges_mod
      USE voidage_mod
      USE deposition_mod
      implicit none
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      integer:: j,i,itmp,narg,len,status
      real(wp),allocatable:: coprxyz(:,:)
      real:: t0,t1
      integer:: ntrial
      character:: c6*6,c5*5
      character(len=132) infile,ctmp
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', ctmp(1:len)
      if(narg < 3)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' natom atom_file ntrial'
         stop
      end if
!
      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) natom
      print *,'natom = ',natom
!
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) ntrial
      print *,'ntrial = ',ntrial
!
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      natom_max = natom
      write (*,*) natom_max
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(atom(1:natom_max))
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
! Read in coordinates and connectivity from xyz file
      read(14,*) itmp;if(itmp /= natom) stop 'itmp /= natom'
      read(14,*)
      do i=1,natom
         read (14,*) ctmp,(rxyz(i,j),j=1,3)
         atom(i) = name2atom(trim(ctmp))
      end do
! Read in coordinates and connectivity from pdb file
!      do i=1,natom
!         read (14,'(a6)',advance='no') c6
!         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j=1,3)
!         atom(i) = name2atom(trim(ctmp))
!      end do
!      proximity = 0
!      do i=1,natom
!         read(14,'(a32)')ctmp
!         do j = 1,ncmax(atom(i))
!            c5 = ctmp(6+5*(j)+1:6+5*(j)+5)
!            read( unit=c5,fmt=* ) proximity(i,j)
!         end do
!      end do
      close(14)
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
!      call assign_charge(1,natom)

!write(*,*) 'sum charge = ',sum(charge(1:natom))
!write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
write(*,*) 'n Si = ',count(atom==iSilicon)
write(*,*) 'n O  = ',count(atom==iOxygen)
write(*,*) 'n OH = ',count(atom==iOxygenH)
write(*,*) 'ntot = ',count(atom==iSilicon) + count(atom==iOxygen) + count(atom==iOxygenH)
!
      call init_probe_mols()
      call INIT_NLIST(boxl,boxl,boxl,5.0_wp/angstrom)
      call NEW_NLIST
!
!call cpu_time(t0)
!      open(unit=30,file='void_points_N2.out')
!      call voidage_profile(probe_N2,ntrial)
!      open(unit=20,file='voidage_N2.out')
!      call write_file(20,voidage)
!      write(*,*) 'N2 voidage = ',sum(voidage)/size(voidage)
!call cpu_time(t1)
!      close(30)
!print *,'time taken for N2 = ',t1-t0
!!
!call cpu_time(t0)
!      open(unit=30,file='void_points_O2.out')
!      call voidage_profile(probe_O2,ntrial)
!      open(unit=21,file='voidage_O2.out')
!      call write_file(21,voidage)
!      write(*,*) 'O2 voidage = ',sum(voidage)/size(voidage)
!call cpu_time(t1)
!      close(30)
!print *,'time taken for O2 = ',t1-t0
!
call cpu_time(t0)
      open(unit=30,file='void_points_He.out')
      call voidage_profile(probe_He,ntrial)
      open(unit=21,file='voidage_He.out')
      call write_file(21,voidage)
      write(*,*) 'He voidage = ',sum(voidage)/size(voidage)
call cpu_time(t1)
      close(30)
print *,'time taken for O2 = ',t1-t0

call cpu_time(t0)
      open(unit=30,file='void_points_CO2.out')
      call voidage_profile(probe_CO2,ntrial)
      open(unit=22,file='voidage_CO2.out')
      call write_file(22,voidage)
      write(*,*) 'CO2 voidage = ',sum(voidage)/size(voidage)
call cpu_time(t1)
      close(30)
print *,'time taken for CO2 = ',t1-t0

contains

   subroutine write_file(iu,voidage)
      integer,intent(in):: iu
      real(wp),intent(in):: voidage(:,:,:)
      integer:: k,j,i
      real(wp):: x,y,z
      do k = 1,ncelz
         z = (k-0.5_wp)*delz - boxl2
         do j = 1,ncely
            y = (j-0.5_wp)*dely - boxl2
            do i = 1,ncelx
               x = (i-0.5_wp)*delx - boxl2
               write(iu,'(3f10.4,g16.8)') x,y,z,voidage(i,j,k)
            end do
         end do
      end do
      write(iu,'(/)')
   end subroutine


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
      INTEGER :: IARGC
      EXTERNAL   IARGC
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

end program simbox


