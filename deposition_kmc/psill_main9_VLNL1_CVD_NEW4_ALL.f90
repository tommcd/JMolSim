
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
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, kboltzT
   integer:: natom,natom_max
   integer:: kmc_step, nrelax
END MODULE global_vars_mod

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
   integer(i4b),public,save:: irand_seed = -1
!
   contains

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

   function rand()
!-----random number generator
      real(wp):: rand
      integer(i4b),parameter::ia=16807,im=2147483647,iq=127773,ir=2836
      integer(i4b),save:: k
      real(wp),save::am
      if (irand_seed <= 0 .or. iy < 0 ) then
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
      if (iy < 0)iy = iy+im
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

!!>include 'files_mod.f90'

MODULE files_mod
   public get_free_file_unit, myflush
CONTAINS

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
   USE global_vars_mod, only: angstrom, K_ev, qstar, pi
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
   integer,parameter:: ntyplj = 9  ! number of LJ types
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
                        eps_O_O2 /)
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
                        sig_O_O2 /)
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
   USE global_vars_mod, only: angstrom,pi
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
! united atom SiOH4
   real(wp),parameter:: r_tet = bondl_SiO + sigma_OH*0.5_wp
!
   real(wp),parameter:: sigLJ(0:3) = (/ sig_O_sil, sig_Si_sil,sigma_OH, sig_H_OH_sil /)
   real(wp),parameter:: sigLJ_2(0:3) = sigma*0.5_wp
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
   integer,parameter:: ncmax(-1:ntyp)=(/0,2,4,2,1,2,1/)
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

END MODULE

!!>include 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:)
!
   integer,parameter:: NDIM = 3          ! Number of dimensions (usually 3)
!
! The periodicity NPER is defined here:
! NPER = 0,  not periodic
! NPER = 1,  periodic in x direction only
! NPER = 2,  periodic in x and y (but NOT z)
! NPER = 3,  periodic in all 3 dimensions
! and box_periodic(i) is set to .TRUE. if dimension i is periodic
!
   integer,private:: i
   integer,parameter:: NPER = 2
   logical,parameter:: box_periodic(ndim) = (/ ( (i <= nper), i=1,ndim ) /)
   real(wp):: boxl(NDIM),boxli(NDIM),boxl2(NDIM)
!
   interface pbc
      module procedure  pbc_v, pbc_a, pbc_a2 , pbc_n_v, pbc_n_a, pbc_n_a2
   end interface
CONTAINS

!   PURE SUBROUTINE pbc_v(r)
!      real(wp),intent(inout):: r(:)
!      r(1:NDIM) = r(1:NDIM) - boxl*anint(r(1:NDIM)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_a(r)
!      real(wp),intent(inout):: r(:,:)
!      r(:,1:NDIM) = r(:,1:NDIM) - boxl*anint(r(:,1:NDIM)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_a2(r)
!      real(wp),intent(inout):: r(:,:,:)
!      r(:,:,1:NDIM) = r(:,:,1:NDIM) - boxl*anint(r(:,:,1:NDIM)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_n_v(r,n)
!      integer,intent(in):: n
!      real(wp),intent(inout):: r(:)
!      r(1:n) = r(1:n) - boxl*anint(r(1:n)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_n_a(r,n)
!      integer,intent(in):: n
!      real(wp),intent(inout):: r(:,:)
!      r(:,1:n) = r(:,1:n) - boxl*anint(r(:,1:n)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_n_a2(r,n)
!      integer,intent(in):: n
!      real(wp),intent(inout):: r(:,:,:)
!      r(:,:,1:n) = r(:,:,1:n) - boxl*anint(r(:,:,1:n)*boxli)
!   END SUBROUTINE

!
   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      where (abs(r(1:NPER)) > boxl2(1:NPER)) &
         r(1:NPER) = r(1:NPER) - sign(boxl(1:NPER),r(1:NPER))
!     if (abs(r(1)) > boxl2(1)) r(1) = r(1) - sign(boxl(1),r(1))
!     if (abs(r(2)) > boxl2(2)) r(2) = r(2) - sign(boxl(2),r(2))
!     if (abs(r(3)) > boxl2(3)) r(3) = r(3) - sign(boxl(3),r(3))
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      integer:: i
      do i = 1,size(r,1)
         where (abs(r(i,1:NPER)) > boxl2(1:NPER)) &
            r(i,1:NPER) = r(i,1:NPER) - sign(boxl(1:NPER),r(i,1:NPER))
!        if (abs(r(i,1)) > boxl2(1)) r(i,1) = r(i,1) - sign(boxl(1),r(i,1))
!        if (abs(r(i,2)) > boxl2(2)) r(i,2) = r(i,2) - sign(boxl(2),r(i,2))
!        if (abs(r(i,3)) > boxl2(3)) r(i,3) = r(i,3) - sign(boxl(3),r(i,3))
      end do
   END SUBROUTINE

   PURE SUBROUTINE pbc_a2(r)
      real(wp),intent(inout):: r(:,:,:)
      integer:: i,j
      do j = 1,size(r,2)
      do i = 1,size(r,1)
         where (abs(r(i,j,1:NPER)) > boxl2(1:NPER)) &
            r(i,j,1:NPER) = r(i,j,1:NPER) - sign(boxl(1:NPER),r(i,j,1:NPER))
!        if (abs(r(i,j,1)) > boxl2(1)) r(i,j,1) = r(i,j,1) - sign(boxl(1),r(i,j,1))
!        if (abs(r(i,j,2)) > boxl2(2)) r(i,j,2) = r(i,j,2) - sign(boxl(2),r(i,j,2))
!        if (abs(r(i,j,3)) > boxl2(3)) r(i,j,3) = r(i,j,3) - sign(boxl(3),r(i,j,3))
      end do
      end do
   END SUBROUTINE

   PURE SUBROUTINE pbc_n_v(r,n)
      integer,intent(in):: n
      real(wp),intent(inout):: r(:)
      where (abs(r(1:n)) > boxl2(1:n)) r(1:n) = r(1:n) - sign(boxl(1:n),r(1:n))
   END SUBROUTINE

   PURE SUBROUTINE pbc_n_a(r,n)
      integer,intent(in):: n
      real(wp),intent(inout):: r(:,:)
      integer:: i
      do i = 1,size(r,1)
         where (abs(r(i,1:n)) > boxl2(1:n)) r(i,1:n) = r(i,1:n) - sign(boxl(1:n),r(i,1:n))
      end do
   END SUBROUTINE

   PURE SUBROUTINE pbc_n_a2(r,n)
      integer,intent(in):: n
      real(wp),intent(inout):: r(:,:,:)
      integer:: i,j
      do j = 1,size(r,2)
      do i = 1,size(r,1)
         where (abs(r(i,j,1:n)) > boxl2(1:n)) r(i,j,1:n) = r(i,j,1:n) - sign(boxl(1:n),r(i,j,1:n))
      end do
      end do
   END SUBROUTINE

   SUBROUTINE WRITE_XYZ(iu)
      USE atom_types_mod, only: atom_name,atom
      integer,intent(in):: iu
      integer:: i,j
      write(iu,*) natom
      write(iu,'(6f14.8)') boxl*Angstrom
      do i = 1,natom
         write(iu,'(a2,6(1x,f14.8))') atom_name(atom(i)),(rxyz(i,j)*Angstrom,j = 1,NDIM)
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu)
      integer,intent(in):: iu
      integer:: i
      real(wp):: r(NDIM)
      character(32):: ctmp
      read(iu,*) natom
      read(iu,*) ctmp
      do i = 1,natom
         read(iu,*) ctmp,r
         rxyz(i,1:NDIM) = r/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod


!!>include 'connectivity_mod.f90'

MODULE connectivity_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom
   USE atom_types_mod
   implicit none
   integer,allocatable:: proximity(:,:)
!
CONTAINS

   SUBROUTINE write_proximity(iu,kmc_step)
!     write topology/connectivity information to file unit iu
      integer,intent(in):: iu,kmc_step
      integer:: i,j,k
      write(iu,*)'# kmc_step = ',kmc_step
      write(iu,*) natom
      do i = 1,natom
         k = count(proximity(i,:) > 0)
         write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
            i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE write_proximity


   SUBROUTINE read_proximity(iu,kmc_step)
!     Read topology/connectivity information from file unit iu
      integer,intent(in):: iu,kmc_step
      integer:: i,j,ic,ii,it
      character(32):: ctmp
      read(iu,*) ctmp,ctmp,ctmp,it
      if (it /= kmc_step)stop 'it /= kmc_step'
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


   PURE FUNCTION is_2Si_ring(iO1)
!    This function will check whether the input oxygen is
!    involved in a 2 Si ring
!    takes a bridging oxygen and checks if the proximity
!    of atom iSi2 contains any atom also bonded to atom
!    iSi1 resulting in a 2Si_ring. i.e. fig B below
!
!   Fig A.   iO1     iO3      Fig B.    iO1     iO3
!            / \    /                   / \    /
!        iSi1   iSi2 - iO2          iSi1   iSi2
!                   \                   \ /    \
!                    iO4                iO2     iO4
!
      integer,intent(in):: iO1
      integer:: iSi1,iSi2,iat,i
      logical:: is_2Si_ring
!
      is_2Si_ring = .FALSE.
      iSi1 = proximity(iO1,1)
      iSi2 = proximity(iO1,2)
      do i = 1, ncmax(atom(iSi2))
         if (proximity(iSi2,i)==iO1) cycle
         iat = proximity(iSi2,i)
         if (ANY(proximity(iat,:)==iSi1)) then
            is_2Si_ring = .TRUE.
            RETURN
         end if
      end do
   END FUNCTION


   PURE FUNCTION OH_groups_OK(iOH1,iOH2)
      logical:: OH_groups_OK
      integer,intent(in):: iOH1,iOH2
      integer:: iSi1,iSi2
      OH_groups_OK = .FALSE.
      iSi1 = proximity(iOH1,1)
      if(atom(iSi1) /= iSilicon) iSi1 = proximity(iOH1,2)
      iSi2 = proximity(iOH2,1)
      if(atom(iSi2) /= iSilicon) iSi2 = proximity(iOH2,2)
      if(iSi1 == iSi2) RETURN
      if (.NOT.nearest_neighbor2(iSi1,iSi2)) OH_groups_OK = .TRUE.
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

   SUBROUTINE check_proximity2(ifirst,ilast,consistent,ierr)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      integer,intent(in):: ifirst,ilast
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k,ia,a1,a2,i1,i2
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         ia = atom(i)
         select case(ia)
         case(iSilicon)  ! Check Si is bonded to 4 atoms
            do j = 1,ncmax(ia)
               k = proximity(i,j)
               if (k == 0) GOTO 100
               if (atom(k) == iSilicon) GOTO 100   ! but not Si
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = proximity(i,1)
            if (i1 == 0) GOTO 100
            a1 = atom(proximity(i,1))
            if (a1 /= iSilicon) GOTO 100
            i2 = proximity(i,2)  !  must be bonded to  2 Silicons
            if (i2 == 0) GOTO 100
            a2 = atom(proximity(i,2))
            if (a2 /= iSilicon) GOTO 100
            if (any(proximity(i,3:) /= 0)) GOTO 100
         case(iOxygenH)  ! OH Oxygen must be bonded to Si
            i1 = proximity(i,1)
            ! if (i1 == 0) GOTO 100
            if (i1 == 0) then
               i1 = proximity(i,2)
               if (i1 == 0) GOTO 100
            end if
            a1 = atom(i1)
            if (a1 /= iSilicon) GOTO 100
            ! if ( any(proximity(i,2:) /= 0) ) GOTO 100
            if ( any(proximity(i,3:) /= 0) ) GOTO 100
         case default
            stop 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ierr = i
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
            write(*,*)'delete_group0: error'
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
            write(*,*)'delete_group2: error'
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

!!>include 'list_mod.f90'

MODULE list_mod
!
   integer,private,parameter:: n_atoml_max = 4000
   TYPE list
      integer:: ind(n_atoml_max) = 0
      integer:: n = 0
      integer:: ip = 0
   END TYPE
!
CONTAINS

   SUBROUTINE add_to_list(lst,iat)
      type(list),intent(inout):: lst
      integer,intent(in):: iat
      lst%n = lst%n + 1
      lst%ind(lst%n) = iat
   END SUBROUTINE

   SUBROUTINE remove_from_list(lst,iat)
      type(list),intent(inout):: lst
      integer,intent(in):: iat
      integer:: i
      do i = 1,lst%n
         if(lst%ind(i)==iat)exit
      end do
      lst%ind(i) = lst%ind(lst%n)
      lst%ind(lst%n) = 0
      lst%n = lst%n - 1
      lst%ip = 0
   END SUBROUTINE

   SUBROUTINE remove_posn_list(lst,ipos)
      type(list),intent(inout):: lst
      integer,intent(in):: ipos
      lst%ind(ipos) = lst%ind(lst%n)
      lst%ind(lst%n) = 0
      lst%n = lst%n - 1
      lst%ip = 0
   END SUBROUTINE

   SUBROUTINE rand_from_list(lst,iat)
      USE rand_mod
      integer,intent(out):: iat
      type(list),intent(inout):: lst
      integer:: j
      if(lst%n > 0)then
         j = int(rand()*lst%n) + 1
         iat = lst%ind(j)
         lst%ip = j
      else
         iat = 0
         lst%ip = 0
      end if
   END SUBROUTINE

END MODULE list_mod

!!>include 'charges_mod.f90'

MODULE charges_mod
   USE precision_mod
   USE atom_types_mod
   USE seaton_mod, only: qi => q_seaton
   implicit none
   PRIVATE
   real(wp),allocatable:: charge(:)
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
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      real(wp):: qo,qh,qoh
      integer:: i
      qoh = qi(iOxygenH) + qi(iHydrogen)
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
   END SUBROUTINE assign_charge

END MODULE charges_mod



!!>include 'bond_list_mod.f90'

MODULE bond_list_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer:: nbondtot
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars_mod, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types_mod, only: atom,ncmax
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
      write(unit = ctmp,fmt = '(i4.4)')iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.pdb'
      inquire(file= trim(frame_file),iostat= ios,exist= lexists,opened = lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit= iu,file= trim(adjustl(frame_file)),form= 'formatted')
   END SUBROUTINE new_frame_file

   SUBROUTINE write_frame(iu,ifirst,ilast)
      USE files_mod
      USE bond_list_mod
      USE global_vars_mod, only: natom,angstrom
      USE connectivity_mod
      USE atom_types_mod, only: atom,atom_name,ncmax
      USE coordinates_mod, only: rxyz,boxl2
      integer,intent(in):: iu,ifirst,ilast
      integer:: i,j,ib,ia,iu2
      logical:: lopen,fix_bonds
      character(len= 132):: frame_file
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
      inquire(unit= iu,name= frame_file)
      j = len_trim(frame_file)
      frame_file = frame_file(1:j-3)//'scr'
      fix_bonds=.false.
      do i = 1,nbondtot
         ia = ibond(1,i)
         ib = ibond(2,i)
         if (ANY( ABS(rxyz(ia,1:3)-rxyz(ib,1:3)) > 10.0_wp/angstrom )) then
            fix_bonds = .true.
            inquire(file= frame_file,opened= lopen)
            if (.not.lopen) then
               call get_free_file_unit(iu2)
               open(unit= iu2,file= trim(frame_file))
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
      write(iu,110)'HETATM',99991,'AR','   ',99991,-boxl2(1)*angstrom,-boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99992,'AR','   ',99992,-boxl2(1)*angstrom, boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2(1)*angstrom,-boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2(1)*angstrom, boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99995,'AR','   ',99995,-boxl2(1)*angstrom,-boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99996,'AR','   ',99996,-boxl2(1)*angstrom, boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2(1)*angstrom,-boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2(1)*angstrom, boxl2(2)*angstrom,5.0
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


!!>include 'nlist_mod2.f90'

MODULE NLIST_MOD
   USE precision_mod, only: wp
   USE coordinates_mod
   implicit none
   PRIVATE
   real(wp),private:: delx,dely,delz,delmin
   real(wp),private:: delxi,delyi,delzi
   integer,public:: ncelx,ncely,ncelz,ncelt
   integer,parameter:: neighmx = 350
   integer,allocatable,public:: hoc(:),hoc_old(:)
!   integer,public:: ncell(neighmx)
   integer,allocatable,public:: ll_old(:),ll(:),lr(:),lr_old(:)
   public:: NEIGCELL,cell,INIT_NLIST,NEW_NLIST,print_cell
   public:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll,store_nlist,restore_nlist
!
!  ll(i)      : linked list particle i
!  hoc(ic)    : head of chain cell ic
!  ncelx,y,z  : number of cells in x, y or z direction
!  ncelt      : total number of cells
!
CONTAINS

   PURE FUNCTION nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   END FUNCTION

   SUBROUTINE INIT_NLIST(Lx,Ly,Lz,Rc,natom_max)
      real(wp),intent(in):: Lx,Ly,Lz,Rc
      integer,intent(in):: natom_max
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max),lr_old(natom_max))
      ncelx = INT(Lx/Rc); ncely = INT(Ly/Rc); ncelz = INT(Lz/Rc)
      ncelt = ncelx*ncely*ncelz
      delx = Lx/ncelx; dely = Ly/ncely; delz = Lz/ncelz
      delmin = min(delx,dely,delz)
      delxi = 1.0_wp/delx; delyi = 1.0_wp/dely; delzi = 1.0_wp/delz
      write(*,'(a,3i6)')'ncelx,y,z = ',ncelx,ncely,ncelz
      write(*,*)'ncelt = ',ncelt
      write(*,'(a,3f12.6)')'delx,y,z = ',delx,dely,delz
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
      CELL = INT((XYZ(1) + boxl2(1))*delxi) &
           + INT((XYZ(2) + boxl2(2))*delyi)*ncelx &
           + INT((XYZ(3) + boxl2(3))*delzi)*ncelx*ncely + 1
   END FUNCTION CELL

   SUBROUTINE NEW_NLIST(ifirst,ilast)
!     makes a new neighbour list using the linked-list algorithm
      integer,intent(in):: ifirst,ilast
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head-of-chain
      ! make linked list:
      do i = ifirst,ilast
         ! determine cell number
         ic = CELL(rxyz(i,:))
!print *,i,rxyz(i,:)
         if (ic < 1 .or. ic > ncelt) stop 'NEW_NLIST: cell number out of bounds'
         ! update linked - list and head of chain
         LL(i) = HOC(ic)
         if (HOC(ic) /= 0) LR(HOC(ic)) = i
         HOC(ic) = i
      end do
   END SUBROUTINE NEW_NLIST

   SUBROUTINE STORE_NLIST
      ll = ll_old
      lr = lr_old
      hoc = hoc_old
   END SUBROUTINE STORE_NLIST

   SUBROUTINE RESTORE_NLIST
      ll_old = ll
      lr_old = lr
      hoc_old = hoc
   END SUBROUTINE RESTORE_NLIST

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
      integer,intent(in):: ic,nlayer
      integer,intent(out):: neigh
      integer,intent(out):: ncell(:)
      integer ix,iy,iz,icx,icy,icz,iccx,iccy,iccz,n2
!
! The periodicity (box_periodic) is set in coordinates_mod
!
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
            if (.NOT.box_periodic(3)) CYCLE
            iccz = iccz + ncelz
         else if (iccz > ncelz) then
            if (.NOT.box_periodic(3)) CYCLE
            iccz = iccz - ncelz
         end if
         do iy = -nlayer,nlayer
            iccy = icy + iy
            if (iccy < 1) then
               if (.NOT.box_periodic(2)) CYCLE
               iccy = iccy + ncely
            else if (iccy > ncely) then
               if (.NOT.box_periodic(2)) CYCLE
               iccy = iccy - ncely
            end if
            do ix = -nlayer,nlayer
               iccx = icx + ix
               if (iccx < 1) then
                  if (.NOT.box_periodic(1)) CYCLE
                  iccx = iccx + ncelx
               else if (iccx > ncelx) then
                  if (.NOT.box_periodic(1)) CYCLE
                  iccx = iccx - ncelx
               end if
               neigh = neigh + 1
               ncell(neigh) = iccx + (iccy - 1)*ncelx + (iccz - 1)*n2
            end do
         end do
      end do
   END SUBROUTINE NEIGCELL

!  SUBROUTINE NEIGCELL_B(ic)
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
      write(iu,*) 'icx,y,z = ',icx,icy,icz
   END SUBROUTINE

   PURE SUBROUTINE Z_NEIGH_CELL(iz,neigh,ncell)
!     determines the neigh cells in the iz'th z-layer
      integer,intent(in):: iz
      integer,intent(out):: neigh
      integer,intent(out):: ncell(:)
      integer:: i
      neigh = ncelx*ncely
      do i = 1,neigh
         ncell(i) = neigh*(iz - 1) + i
      end do
   END SUBROUTINE

   SUBROUTINE ixyzcell(ic,ix,iy,iz)
      integer,intent(in):: ic
      integer,intent(out):: ix,iy,iz
      integer:: n2
      n2 = ncelx*ncely
      iz = ic/n2 + 1
      if (mod(ic,n2) == 0) iz = iz - 1
      iy = (ic - (iz - 1)*n2)/ncelx + 1
      if (mod(ic - (iz - 1)*n2,ncelx) == 0) iy = iy - 1
      ix = ic - (iy - 1)*ncelx - (iz - 1)*n2
   END SUBROUTINE ixyzcell

   PURE FUNCTION icell(ix,iy,iz)
      integer:: icell
      integer,intent(in):: ix,iy,iz
      icell = (iz - 1)*ncelx*ncely + (iy - 1)*ncelx + ix
   END FUNCTION icell

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

END MODULE NLIST_MOD


!!>include 'verlet_list_nl.f90'

MODULE verlet_list_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom_max,natom
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nabormx = 100
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
      integer,parameter:: nl_max = 2
      integer:: i,j,ic,jj,nc,nl,nn,neighbors((2*nl_max+1)**3)
      nl = nlayers_ll(rv - tiny(1.0))
      if(nl > nl_max) STOP ' NEW_VLIST: nl > nl_max'
      nlist(1:natom)=0
      rxyz_old(1:natom,1:NDIM) = rxyz(1:natom,1:NDIM)
      call NEW_NLIST(1,natom)
      do i = 1, natom
      ic = CELL(rxyz(i,1:NDIM))  ! link list
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
         nc = neighbors(jj)
         if (nc == 0)cycle cell_loop
         j = hoc(nc)
         do while(j /= 0)
            if (j /= i .AND. j > i) then
            r = rxyz(j,:) - rxyz(i,:)
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
      if (ANY(NLIST(1:natom) >= nabormx) )stop 'NLIST(:) >= nabormx'
!!do i = 1,natom
!!   call qsort(nlist(i),list(i,1:nlist(i)),cclst)
!!   list(i,1:nlist(i))=cclst(1:nlist(i))
!!end do
   END SUBROUTINE NEW_VLIST

   logical function update(i)
      integer,intent(in):: i
      real(wp):: r(3)
      r(1:NDIM) = rxyz(i,1:NDIM)-rxyz_old(i,1:NDIM)
      call pbc(r)
      update = (dot_product(r,r) > SKIN2*0.25_wp)
   end function

   FUNCTION UPDATE2(ifirst,ilast)
!     decides whether the list needs to be reconstructed
      logical:: UPDATE2
      integer,intent(in):: ifirst,ilast
      real(wp):: dispmx
!     a conservative test of the list skin crossing
      dispmx = MAXVAL(ABS( rxyz(ifirst:ilast,1:NDIM) - rxyz_old(ifirst:ilast,1:NDIM) ))
      dispmx = 2.0_wp*SQRT( 3.0_wp*dispmx**2 )
      UPDATE2 = ( dispmx > skin )
   END FUNCTION

END MODULE verlet_list_mod


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
   real(wp),parameter:: kbond(1:3)=(/ KSiO,KSiSi,KSiO /)
   real(wp),parameter:: abond(1:3)=(/ ASiO,ASiSi,ASiO /)
   real(wp),parameter:: acosang(0:2)=(/ 1.0_wp,1.0_wp/3.0_wp,1.0_wp /)
END MODULE Keating_parameters


!!>include 'repul_force_NL_vonAlfthan.f90'

MODULE repul_energy_mod
   USE precision_mod, only: wp
   USE global_vars_mod
   implicit none
   real(wp),parameter,private:: rden = 2.6_wp/Angstrom
   real(wp),parameter,private:: rden2 = rden**2
   real(wp),parameter,private:: ff = 0.8_wp*(Angstrom**4)   ! 8000.0_wp eV nm^-4
   real(wp):: repulenergy
CONTAINS

   SUBROUTINE repul_force(ifirst,ilast,en)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp),intent(out):: en
      real(wp):: rr(NDIM),r2,r(NDIM)
      integer:: M,L,ic,nc,inn,nn,neigbors(27)
!
      en = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      r = rxyz(L,:)
      ic = CELL(r) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            rr(1:NDIM) = rxyz(M,1:NDIM) - r(1:NDIM)
            call pbc(rr)
            r2 = dot_product(rr,rr)
            if (r2 <= rden2) then
               if (nearest_neighbor2(L,M)) GOTO 100
               en = en + (r2 - rden2)**2
               fxyz(L,:) = fxyz(L,:) + 2.0_wp*ff*(r2 - rden2)*(rr)
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      en = en*0.5_wp*ff
   END SUBROUTINE repul_force

   SUBROUTINE repul_force2(ifirst,ilast,en)
      USE coordinates_mod
      USE atom_types_mod, only: atom
      USE connectivity_mod
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp),intent(out):: en
      real(wp):: r2,dx,dy,dz,x,y,z
      integer:: M,L,ic,nc,inn,nn,neigbors(27),i,j,k
!
      en = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      !r = rxyz(L,:)
      x = rxyz(L,1)
      y = rxyz(L,2)
      z = rxyz(L,3)
      ic = CELL((/x,y,z/)) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            dx = rxyz(M,1) - x
            dy = rxyz(M,2) - y
            dz = rxyz(M,3) - z
            if (abs(dx) > boxl2(1)) dx = dx - sign(boxl(1),dx)
            if (abs(dy) > boxl2(2)) dy = dy - sign(boxl(2),dy)
            if (abs(dz) > boxl2(3)) dz = dz - sign(boxl(3),dz)
            r2 = dx*dx + dy*dy + dz*dz
            if (r2 < rden2) then
               do i = 1,ncmax(atom(L))
                  K = proximity(L,i)
                  if (K == 0) cycle
                  if (K == M) GOTO 100
                  do j = 1,ncmax(atom(K))
                     if (proximity(K,j) == M) GOTO 100
                  end do
               end do
               en = en + 0.5_wp*ff*(r2 - rden2)**2
               fxyz(L,1) = fxyz(L,1) + 2.0_wp*ff*(r2 - rden2)*dx
               fxyz(L,2) = fxyz(L,2) + 2.0_wp*ff*(r2 - rden2)*dy
               fxyz(L,3) = fxyz(L,3) + 2.0_wp*ff*(r2 - rden2)*dz
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
   END SUBROUTINE repul_force2


   PURE FUNCTION repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom
      USE connectivity_mod
      USE NLIST_MOD
      real(wp):: repul_energy
      integer,intent(in):: ifirst,ilast
      real(wp):: r2,dx,dy,dz,x,y,z,en
      integer:: M,L,ic,nc,inn,nn,neigbors(27),i,j,k
!
      en = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      x = rxyz(L,1)
      y = rxyz(L,2)
      z = rxyz(L,3)
      ic = CELL((/x,y,z/))
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            dx = rxyz(M,1) - x
            dy = rxyz(M,2) - y
            dz = rxyz(M,3) - z
            if (abs(dx) > boxl2(1)) dx = dx - sign(boxl(1),dx)
            if (abs(dy) > boxl2(2)) dy = dy - sign(boxl(2),dy)
            if (abs(dz) > boxl2(3)) dz = dz - sign(boxl(3),dz)
            r2 = dx*dx + dy*dy + dz*dz
            if (r2 < rden2) then
               do i = 1,ncmax(atom(L))
                  K = proximity(L,i)
                  if (K == 0) cycle
                  if (K == M) GOTO 100
                  do j = 1,ncmax(atom(K))
                     if (proximity(K,j) == M) GOTO 100
                  end do
               end do
               en = en + (r2 - rden2)**2
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      repul_energy = 0.5_wp*ff*en
   END FUNCTION repul_energy

   PURE FUNCTION repul_energy_NL(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      real(wp):: repul_energy_NL
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2,en
      integer:: M,L,ic,nc,inn,nn,neigbors(125)
      en = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      ic = CELL(rxyz(L,:)) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            rr(1:NDIM) = rxyz(M,1:NDIM) - rxyz(L,1:NDIM)
            call pbc(rr)
            r2 = dot_product(rr,rr)
            if (r2 <= rden2)then
               en = en + (r2 - rden2)**2
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      repul_energy_NL = en*0.5_wp*ff
   END FUNCTION repul_energy_NL

   PURE FUNCTION repul_energy_VL(ifirst,ilast)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE verlet_list_mod
      real(wp):: repul_energy_VL
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2,en
      integer:: M,L,jj
      en = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
         atom_loop: do jj = 1, NLIST(L)  ! Verlet-list
            M = LIST(L,jj)
            if (nearest_neighbor2(L,M)) CYCLE atom_loop
            if (M == L) CYCLE atom_loop
            rr(1:NDIM) = rxyz(M,1:NDIM) - rxyz(L,1:NDIM)
            call pbc(rr)
            r2 = rr(1)**2+rr(2)**2+rr(3)**2
            if (r2 <= rden2) then
               en = en + (r2 - rden2)**2
            end if
         end do atom_loop
      end do outer_atom_loop
      repul_energy_VL = en*0.5_wp*ff
   END FUNCTION repul_energy_VL

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
      USE atom_types_mod, only: atom,ncmax,iHydrogen
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
      do i = ifirst,ilast
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

!!>include 'HKNonLattice2.f90'

MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster=0,n_cluster_old=0
   integer,parameter,private :: ncmx = 2000
   integer,allocatable:: atomL(:),atomL_old(:)
   integer,allocatable,private:: cluster(:,:)
   integer,private:: clusterC(ncmx)
CONTAINS

   SUBROUTINE Init_HKNonLattice(natom_max)
      integer,intent(in):: natom_max
      if(.not.allocated(atomL))then
         allocate(atomL(natom_max),atomL_old(natom_max))
         allocate(cluster(ncmx,natom_max))
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

! some routines used mainly for debugging
   SUBROUTINE analyse_cluster(natom)
      implicit none
      integer,intent(in):: natom
      integer:: i,L
      clusterC(1:n_cluster)=0
      cluster(1:n_cluster,1:ncmx)=0
      do i = 1,natom
         L = AtomL(i)
         clusterC(L) = clusterC(L) + 1
         cluster(L,clusterC(L))=i
      end do
      ! some consistency checks
      do i = 1,n_cluster
         if (clusterC(i) /= count(atomL(1:natom)==i)) then
            stop 'error in analyse_cluster'
         end if
      end do
      if (sum(clusterC(1:n_cluster))/=natom) then
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
         if (clusterC(i) /= count(atomL(1:natom)==i)) then
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
   type(probe_mol),save:: probe_SiOH4, probe_SiO4, probe_tip3p
   type(probe_mol),save:: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet
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
      if(present(radius)) P%rad = radius
      if(present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types_mod
      USE seaton_mod
      USE charges_mod
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigma_H2O*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom = (/ iO_H2O, iH_H2O, iH_H2O /)
      probe_TIP3P%rad = sigi2(probe_TIP3P%atom)
!      probe_TIP3P%q = q_seaton(probe_TIP3P%atom)
      probe_TIP3P%q = (/q_O_H2O, -0.5_wp*q_O_H2O, -0.5_wp*q_O_H2O /)
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
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigi2(probe_SiO4%atom)
! O2
      call new_probe_mol(probe_O2,2,sig_O_O2*0.5_wp)
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_O2*0.5_wp /)
! N2
      call new_probe_mol(probe_N2,2,sig_N_N2*0.5_wp)
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_N2*0.5_wp /)
! CO2
      call new_probe_mol(probe_CO2,3,sig_O_CO2*0.5_wp)
      probe_CO2%rad(1) = sig_C_CO2*0.5_wp
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp,-bondl_CO2 /)
! Si(OH)4
      write(*,*) 'SiOH4'
      write(*,*) probe_SiOH4%n
      write(*,'(9i6)') probe_SiOH4%atom
      write(*,'(9f12.6)') probe_SiOH4%rad
      write(*,'(9f12.6)') probe_SiOH4%q
      write(*,'(3f12.6)') probe_SiOH4%r
! H2O
      write(*,*) 'H2O TIP3P'
      write(*,*) probe_TIP3P%n
      write(*,'(9i6)') probe_TIP3P%atom
      write(*,'(9f12.6)') probe_TIP3P%rad
      write(*,'(9f12.6)') probe_TIP3P%q
      write(*,'(3f12.6)') probe_TIP3P%r
! SiO4
      write(*,*) 'SiO4'
      write(*,*) probe_SiO4%n
      write(*,'(9f12.6)') probe_SiO4%rad
      write(*,'(3f12.6)') probe_SiO4%r
! O2
      write(*,*) 'O2'
      write(*,*) probe_O2%n
      write(*,'(9f12.6)') probe_O2%rad
      write(*,'(3f12.6)') probe_O2%r
! N2
      write(*,*) 'N2'
      write(*,*) probe_N2%n
      write(*,'(9f12.6)') probe_N2%rad
      write(*,'(3f12.6)') probe_N2%r
! CO2
      write(*,*) 'CO2'
      write(*,*) probe_CO2%n
      write(*,'(9f12.6)') probe_CO2%rad
      write(*,'(3f12.6)') probe_CO2%r
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
      USE atom_types_mod
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

!!>include 'lj_el_mod.f90'

MODULE lj_el_mod
    USE precision_mod, only: wp
    USE global_vars_mod
    implicit none
    private
    save
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: uljcut(:,:),uljsf(:,:)
    real(wp),parameter:: rcut = 12.0_wp/angstrom
    real(wp),parameter:: rcut2 = rcut**2
    PUBLIC:: LJ_INIT,ENERGY_LJ_EL
CONTAINS

   SUBROUTINE LJ_INIT
      USE seaton_mod, only: epsi,sigi,ntyplj
      real(wp):: eij,rstij,rcuti,rcuti6
      integer:: i,j
      rcuti = (1.0_wp/rcut)
      rcuti6 = rcuti**6
      if(.not.allocated(aij)) allocate( aij(0:ntyplj,0:ntyplj) )
      if(.not.allocated(bij)) allocate( bij(0:ntyplj,0:ntyplj) )
      if(.not.allocated(uljcut)) allocate( uljcut(0:ntyplj,0:ntyplj) )
      if(.not.allocated(uljsf))  allocate(  uljsf(0:ntyplj,0:ntyplj) )
      do j = 0,ntyplj    ! lj aij,bij parameters & shifted force terms
      do i = 0,ntyplj
        eij = sqrt(epsi(i)*epsi(j))
        rstij = (sigi(i)+sigi(j))*0.5
        aij(i,j) = 4.0*(rstij**12)*eij
        bij(i,j) = 4.0*(rstij**6)*eij
        uljcut(i,j) = rcuti6*(aij(i,j)*rcuti6-bij(i,j))
        uljsf(i,j) = -rcuti6*(12.0*aij(i,j)*rcuti6-6.0*bij(i,j))*rcuti
        write(*,'(2i6,4g18.9)')i,j,aij(i,j),bij(i,j),uljcut(i,j),uljsf(i,j)
      end do
      end do
   END SUBROUTINE

   PURE SUBROUTINE ENERGY_LJ_EL(pr,sysfirst,syslast,Ulj,Uel)
      USE coordinates_mod
      USE atom_types_mod, only: atom
      USE charges_mod, only: charge
      USE probe_mol_mod
      type(probe_mol),intent(in):: pr
      integer,intent(in):: sysfirst,syslast
      real(wp),intent(out):: Ulj,Uel
      real(wp):: rr(3),r2,r1,rr6,dele,qL
      integer:: L,M,iaL
      Ulj = 0.0_wp
      Uel = 0.0_wp
      outer_atom_loop: do L = 1,pr%n
      iaL = pr%atom(L)
      qL = pr%q(L)
      atom_loop: do M = sysfirst,syslast
         rr(1:3) = rxyz(M,1:3) - pr%r(1:3,L)
         call pbc(rr)
         !r2 = dot_product(rr,rr)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         if (r2 < rcut2 .and. M > 0) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(iaL,atom(M))*rr6-bij(iaL,atom(M))) !&
                 !- uljcut(iaL,atom(M)) - uljsf(iaL,atom(M))*(r1-rcut)
            Ulj = Ulj + dele
         end if
         Uel = Uel + qL*charge(M)/r1
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE

   PURE SUBROUTINE EnergyLJEL(atomfirst,atomlast,sysfirst,syslast,Uljel)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE atom_types_mod, only: atom
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
         call pbc(rr)
         r2 = rr(1)**2+rr(2)**2+rr(3)**2
         r1 = sqrt(r2)
         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(atom(L),atom(M))*rr6-bij(atom(L),atom(M))) !&
                 !- uljcut(atom(L),atom(M)) - uljsf(atom(L),atom(M))*(r1-rcut)
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
    USE global_vars_mod
    implicit none
    PUBLIC:: Henrys_law_calc_z
CONTAINS

   SUBROUTINE Henrys_law_calc_z(ntrial,zlower,zupper,pr,fover,Uljel,Khenry,facc)
      USE coordinates_mod
      USE atom_types_mod
!     USE constants_mod, only: sigma_2
      USE seaton_mod, only: sigi2
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE global_vars_mod, only: kboltzT
      USE lj_el_mod
      USE probe_mol_mod
      integer,intent(in):: ntrial
      real(wp),intent(in):: zlower,zupper
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
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2(1)
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2(2)
         z = zlower + rand()*(zupper - zlower)
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
         do i = 1,nat
            call pbc(pr1%r(1:3,i))
         end do
! first check if there is a large overlap
         if(overlap(fover)) cycle insert_loop
         rnaccept = rnaccept + 1.0_wp
         call ENERGY_LJ_EL(pr1,1,natom,dUlj,dUel)
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
   END SUBROUTINE Henrys_law_calc_z


   SUBROUTINE Henrys_law_calc(ntrial,xl,xu,yl,yu,zl,zu,pr,fover,Uljel,Khenry,facc)
      USE coordinates_mod
      USE atom_types_mod
!     USE constants_mod, only: sigma_2
      USE seaton_mod, only: sigi2
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE global_vars_mod, only: kboltzT
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
         do i = 1,nat
            call pbc(pr1%r(1:3,i))
         end do
! first check if there is a large overlap
         if(overlap(fover)) cycle insert_loop
         rnaccept = rnaccept + 1.0_wp
         call ENERGY_LJ_EL(pr1,1,natom,dUlj,dUel)
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
         integer:: k,jj,j,ic,nc,ncell(125),neigh,sig2j
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
               !if (atom(j) /= iSilicon) then
               dr(1:3) = pr1%r(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               sig2j = sigi2(atom(j))
               if (atom(j) == iSilicon) sig2j = 0.6_wp/angstrom
               if (dot_product(dr,dr) < ( fr*( pr1%rad(k) + sig2j ) )**2) then
                  overlap = .true.
                  return
               end if
               !end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
         end do probe_atom_loop
      end function overlap
   END SUBROUTINE

END MODULE Henrys_law_calc_mod

!!>include 'voidage4_mod.f90'

MODULE voidage_mod
   USE precision_mod
   implicit none
CONTAINS

   FUNCTION voidage_calc(ntrial,zlower,zupper,pr)
      USE coordinates_mod
      USE atom_types_mod, only: atom,iSilicon
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
      integer:: i,k,j,ic,jj,nc,it,nat,neigh,ncell(125)
      logical:: rotate
      nat = pr%n
      allocate(probe(3,nat))
      rotate = ((size(probe,2) > 1) .and. nat > 1)
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2(1)
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2(2)
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
                  if (dot_product(dr,dr) < (pr%rad(k) + sigLJ_2(atom(j)))**2) then
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
      USE atom_types_mod, only: atom,iSilicon
      USE constants_mod, only: sigma_2
      USE rand_mod
      USE nlist_mod
      real(wp):: voidage_calc1
      integer,intent(in):: ntrial
      real(wp),intent(in):: zlower,zupper,rprobe
      real(wp):: r(3),dr(3),rnaccept
      integer:: j,ic,jj,nc,it,neigh,ncell(125)
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         r(1) = ( 1.0_wp - 2.0_wp*rand() )*boxl2(1)
         r(2) = ( 1.0_wp - 2.0_wp*rand() )*boxl2(2)
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
                  if (dot_product(dr,dr) < (rprobe + sigLJ_2(atom(j)))**2) then
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

!!>include 'density_profile_mod.f90'

MODULE density_profile_mod
   USE precision_mod, only: wp
   USE atom_types_mod
   implicit none
   integer,parameter,private:: nbinmax = 1000
   integer,private:: nbz,nbx,nby
   real(wp):: ndenz(nbinmax,0:ntyp),nden(10,10,nbinmax,0:ntyp)
   real(wp),private:: delz,delx,dely

CONTAINS

   SUBROUTINE den_profile_z(ibegin,iend,dl)
      USE coordinates_mod
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,delzi
      integer:: j,ii,jj
      top_atom = maxval(rxyz(ibegin:iend,3))
      nbz = int(top_atom/dl) + 1
      delz = dl
      delzi = 1.0_wp/delz
      ndenz(1:nbz,:) = 0.0_wp
      do j = ibegin,iend
         ii = int(rxyz(j,3)*delzi + 1.0_wp)
         if (ii > nbz) then
            if ( rxyz(j,3) == top_atom) then
               ii = nbz
            else
               stop 'ii > nbz'
            end if
         end if
         jj = atom(j)
         ndenz(ii,jj) = ndenz(ii,jj) + 1.0_wp
      end do
      ndenz(1:nbz,:) = ndenz(1:nbz,:)/(boxl(1)*boxl(2)*delz*angstrom**3)
   END SUBROUTINE

   SUBROUTINE den_profile(ibegin,iend,dl)
      USE coordinates_mod
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,delzi,delxi,delyi
      integer:: j,iz,jj,ix,iy
      top_atom = maxval(rxyz(ibegin:iend,3))
      nbx = int(boxl(1)/dl)
      nby = int(boxl(2)/dl)
      delx = boxl(1)/nbx
      dely = boxl(2)/nby
      delz = delx
      nbz = int(top_atom/delz) + 1
      delxi = 1.0_wp/delx
      delyi = 1.0_wp/dely
      delzi = 1.0_wp/delz
      nden = 0.0_wp
      do j = ibegin,iend
         ix = int((rxyz(j,1)+boxl2(1))*delxi) + 1
         iy = int((rxyz(j,2)+boxl2(2))*delyi) + 1
         iz = int(rxyz(j,3)*delzi) + 1
         if (iz > nbz) then
            if (rxyz(j,3) == top_atom) then
               iz = nbz
            else
               stop 'iz > nbz'
            end if
         end if
         jj = atom(j)
         nden(ix,iy,iz,jj) = nden(ix,iy,iz,jj) + 1.0_wp
      end do
      nden(1:nbx,1:nby,1:nbz,:) = nden(1:nbx,1:nby,1:nbz,:)/(delx*dely*delz*angstrom**3)
      forall(iz=1:nbz,jj=0:ntyp) ndenz(iz,jj) = sum(nden(:,:,iz,jj))/(nbx*nby)
   END SUBROUTINE

   SUBROUTINE write_den_profile_z(iu,nstep)
      integer,intent(in):: iu,nstep
      integer:: i,j
      write(iu,'(a,i6)') '#nstep ',nstep
      write(iu,'(a)') '#density [atoms/A^3]'
      do i = 1,nbz
         write(iu,'(10f12.6)') (i-0.5_wp)*delz,(ndenz(i,j),j=0,ntyp)
      end do
      write(iu,'(/)')
   END SUBROUTINE

END MODULE density_profile_mod

!!>include 'deposition4_mod.f90'

MODULE deposition_mod
   USE Henrys_law_calc_mod
   USE precision_mod, only: wp
   implicit none
   logical,parameter:: use_voidage_calc = .TRUE.
   integer,parameter,private:: ntrial = 10000
   integer,parameter,private:: nbzmx = 200
   integer,parameter,private:: nbxmx = 15,nbymx = 15
   integer:: nbin,itop_nonp,itop_nonp_SiO4,itop_nonp_H2O
   real(wp):: nden(nbzmx,0:4),voidage(nbzmx)
   real(wp):: henryk(nbzmx,2),inthenryk(nbzmx,2)
   real(wp):: henryk2(nbxmx,nbymx,nbzmx,4)
   !real(wp):: Eljel(nbxmx,nbxmx,nbzmx,4)
   real(wp):: intvoidage(nbzmx),delz
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


   SUBROUTINE henry_profile_z(probe,dl,ip)
      USE voidage_mod
      USE coordinates_mod
      USE probe_mol_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: dl
      integer,intent(in):: ip
      real(wp):: top_atom,zl,zu,facc,uljel,khenry,fover
      integer:: i,ntrial
      top_atom = maxval(rxyz(1:natom,3))
      nbin = int(top_atom/dl) + 1
      delz = dl
      ntrial = 36000
      fover = 0.8_wp
      do i = 1,nbin
         zl = (i-1)*delz
         zu = i*delz
         call Henrys_law_calc_z(ntrial,zl,zu,probe,fover,Uljel,Khenry,facc)
         henryk(i,ip) = Khenry
      end do
   END SUBROUTINE

   SUBROUTINE henry_profile(probe,dl,ip)
      USE voidage_mod
      USE coordinates_mod
      USE probe_mol_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: dl
      integer,intent(in):: ip
      integer,parameter:: ntrial = 1000
      real(wp),parameter:: fover = 0.8_wp
      real(wp):: facc,uljel,khenry
      real(wp):: xl,xu,yl,yu,zl,zu
      real(wp):: delz,delx,dely,top_atom
!     real::t0,t1
      integer:: i,j,k,nbx,nby,nbz
      top_atom = maxval(rxyz(1:natom,3))
      nbx = int(boxl(1)/dl)
      nby = int(boxl(2)/dl)
      !print *,'nbx = ',nbx
      !print *,'nby = ',nby
      delx = boxl(1)/nbx
      dely = boxl(2)/nby
      delz = delx
      nbz = int(top_atom/delz) + 1
      do k = 1,nbz
         zl = (k-1)*delz
         zu = k*delz
         do j = 1,nby
            yl = (j-1)*dely - boxl2(2)
            yu = j*dely - boxl2(2)
            do i = 1,nbx
               xl = (i-1)*delx - boxl2(1)
               xu = i*delx - boxl2(1)
               !call cpu_time(t0)
               call Henrys_law_calc(ntrial,xl,xu,yl,yu,zl,zu,probe,fover,Uljel,Khenry,facc)
               !call cpu_time(t1)
               henryk2(i,j,k,ip) = Khenry
               !Eljel(i,j,k,ip) = Uljel
               !print '(3i5,4g16.8)',k,j,i,Khenry,Uljel,facc,t1-t0
            end do
         end do
      end do
      forall(k=1:nbz) henryk(k,ip) = sum(henryk2(1:nbx,1:nby,k,ip))/(nbx*nby)
      nbin = nbz
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
         zl = (i-1)*delz
         zu = i*delz
         if (probe%n==1) then
            voidage(i) = voidage_calc1(ntrial,zl,zu,probe%rad(1))
         else
            voidage(i) = voidage_calc(ntrial,zl,zu,probe)
         end if
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


   SUBROUTINE write_voidage_dist(iu,itop_nonp)
      integer,intent(in):: iu,itop_nonp
      integer:: i
      real(wp):: rsum
      rsum = sum(voidage(itop_nonp:nbin))
      write(iu,'(a)') '# voidage'
      do i = itop_nonp,nbin
         write(iu,'(3f12.6)') (i-0.5_wp)*delz,voidage(i)/rsum,intvoidage(i)
      end do
      write(iu,'(/)')
   END SUBROUTINE

   SUBROUTINE write_voidage(iu)
      integer,intent(in):: iu
      integer:: i
      write(iu,'(a)') '# voidage'
      do i = 1,nbin
         write(iu,'(2f12.6)') (i-0.5_wp)*delz,voidage(i)
      end do
      write(iu,'(/)')
      call flush(iu)
   END SUBROUTINE

   SUBROUTINE rand_pos(itop_nonp,z_dep)
      USE rand_mod
      integer,intent(in):: itop_nonp
      real(wp),intent(out):: z_dep
      integer:: iv1(1)
      real(wp):: rnd
      rnd = rand()
      iv1 = MINLOC(ABS(intvoidage(itop_nonp:nbin)-rnd))+itop_nonp-1
!      z_dep = (iv1(1)-0.5)*delz
      z_dep = (iv1(1)-1)*delz + rand()*delz
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

!!>include 'rdf_mod.f90'

MODULE RDF_MOD
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nbin=200
   real(wp),private:: RDF(nbin)
   real(wp),private:: DL1
   PUBLIC:: RDF_CALC, RDF_PRINT
CONTAINS

   SUBROUTINE RDF_CALC
      real(wp):: rv(3),r,r2,rmax,rmax2
      integer:: i,j,k
      rmax = maxval(boxl)*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDF = 0.0_wp
      do i = 1,natom-1
         do j = i + 1,natom
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rmax2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1)+1
            RDF(k) = RDF(k) + 1.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC

   SUBROUTINE RDF_PRINT(io,kmc_step)
      integer,intent(in):: io,kmc_step
      integer:: i
      write(io,*)'# kmc_step = ',kmc_step
      write(io,*)'# (natom) = ',natom
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i-0.5_wp)*DL1,RDF(i)/(natom)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT

END MODULE RDF_MOD

!!>include 'relax_VL_mod.f90'

MODULE relax_mod
   USE precision_mod, only: wp
   implicit none
   public:: RELAX, RELAX_LIST
CONTAINS

   SUBROUTINE RELAX(nr,ifirst,ilast,kbT)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE Keating
      USE repul_energy_mod
      USE verlet_list_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,ifirst,ilast
      real(wp),intent(in):: kbT
!     real(wp),parameter:: de = 0.006_wp
      real(wp):: de = 0.08_wp
      integer:: ir,ia,j,k,nfail,ntot
      real(wp):: coprxyz(3),energy1,energy2,del_energy,dU_kT
!
      nfail = 0
      ntot = 0
      DO ir = 1,NR
      DO ia = ifirst,ilast
         if (atom(ia) == iHydrogen) CYCLE
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
         energy1 = energy1 + repul_energy(ia,ia)

         rxyz(ia,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
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
         energy2 = energy2 + repul_energy(ia,ia)
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
      END DO
!write(888, *) energy1,energy2,del_energy
      END DO
      if (real(nfail)/ntot < 0.5_wp) then
         de = de*1.02_wp
      else
         de = de/1.02_wp
      end if
!write(999,'(f0.6,2(i7,1x),f0.6)') de,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX


   SUBROUTINE RELAX_LIST(nr,nl,alist,kbT)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE Keating
      USE repul_energy_mod
      USE verlet_list_mod
      USE rand_mod
      implicit none
      integer,intent(in):: nr,nl,alist(:)
      real(wp),intent(in):: kbT
!     real(wp),parameter:: de = 0.006_wp
      real(wp):: de = 0.08_wp
      integer:: ir,ia,ii,j,k,nfail,ntot
      real(wp):: coprxyz(3),energy1,energy2,del_energy,dU_kT
!
      nfail = 0
      ntot = 0
      DO ir = 1,NR
      DO ii = 1,nl
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
         energy1 = energy1 + repul_energy(ia,ia)

         rxyz(ia,1) = coprxyz(1) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,2) = coprxyz(2) + (2.0_wp*rand() - 1.0_wp)*de
         rxyz(ia,3) = coprxyz(3) + (2.0_wp*rand() - 1.0_wp)*de
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
         energy2 = energy2 + repul_energy(ia,ia)
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
      END DO
!write(888, *) energy1,energy2,del_energy
      END DO
      if (real(nfail)/ntot < 0.5_wp) then
         de = de*1.02_wp
      else
         de = de/1.02_wp
      end if
!write(999,'(f0.6,2(i7,1x),f0.6)') de,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX_LIST

END MODULE relax_mod

!!>include 'select_atom_mod.f90'

MODULE select_atom_mod
   USE precision_mod
   USE global_vars_mod
   USE constants_mod
   USE atom_types_mod
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
      iatom = arraybonds( int(rand()*n)+1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment3(iz,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      integer,intent(out):: iatom
      integer:: arraybonds(natom_max),n,L,ic,nc,neigh,ncell(125)
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            n = n + 1
            arraybonds(n) = L
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      iatom = arraybonds( int(rand()*n)+1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment2(iz,diam,bondl,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      real(wp),intent(in):: diam,bondl
      integer,intent(out):: iatom
      integer,parameter:: ntrial_max = 100
      real(wp):: r3v(3),alfa,beta
      integer:: arraybonds(natom_max),i,n,L,it,ic,nc,neigh,ncell(125)
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            trial_loop: do it = 1,ntrial_max
               alfa = acos(2.0_wp*rand() - 1.0_wp)
               beta = rand()*2.0_wp*pi
               rxyz(natom+1,3) = rxyz(L,3) + bondl*cos(beta)
               rxyz(natom+1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
               rxyz(natom+1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
               call pbc(rxyz(natom+1,:))
               do i = 1,natom
                  if (atom(i) == iHydrogen) cycle
                  if (atom(i) == iSilicon) cycle
                  if (i == L) cycle
                  r3v(1:3) = rxyz(i,1:3) - rxyz(natom+1,1:3)
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
      iatom = arraybonds( int(rand()*n)+1 )
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
            rxyz(natom+1,3) = rxyz(L,3) + bondl*cos(beta)
            rxyz(natom+1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
            rxyz(natom+1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
            call pbc(rxyz(natom+1,:))
            do i = 1,natom
               if (atom(i) == iSilicon) cycle
               if (i == L) cycle
               r3v(1:3) = rxyz(i,1:3) - rxyz(natom+1,1:3)
               call pbc(r3v)
               if (dot_product(r3v,r3v) <= diam**2) CYCLE trial_loop
            end do
            exit trial_loop
         end do trial_loop
         if (it >= ntrial_max)cycle atom_loop
         n = n + 1
         arraybonds(n) = L
      end do atom_loop
      iv = MINLOC(ABS(rxyz(arraybonds(1:n),3)-z_dep))
      iatom = arraybonds(iv(1))
   END SUBROUTINE get_atom_for_attachment

END MODULE select_atom_mod


!!>include 'attachment_NL_mod.f90'

MODULE attachment_mod
   USE precision_mod
   USE global_vars_mod
   USE constants_mod
   USE atom_types_mod
   USE coordinates_mod
   USE connectivity_mod
   USE rand_mod
   implicit none
CONTAINS

   SUBROUTINE attach_tetrahedra(L,beta,bondl,overlap_check,success)
      USE matvec3d_mod,only: len_3d
      USE nlist_mod
      USE HKNonLattice_mod
      USE rotate_axis_mod
      integer,intent(in):: L
      real(wp),intent(in):: beta,bondl
      logical,intent(in):: overlap_check
      logical,intent(out):: success
      integer,parameter:: ntry=100
      real(wp):: xx,yy,zz,phi,cphi,sphi,cbeta,sbeta
      real(wp):: sth,cth,rxy
      real(wp):: r3v(3),r2,r12(3),r23(3)
      real(wp):: tetra(3,3),aa(3,3)
      integer:: j,LSi,ii,jj,LS,M,ic,nc,nl,neigh,ncell(125)
      logical:: overlap
      nl = nlayers_ll(sigLJ(iOxygenH))
      cbeta = cos(beta)
      sbeta = sin(beta)
      success = .false.
! Attaching a new Si atom
      natom = natom + 1
      LS = natom
      LSi = 0
      do j = 1,ncmax(atom(L))
         LSi = proximity(L,j)
         if(LSi /= 0)then
            if( atom(LSi) == iSilicon ) EXIT
         end if
      end do
if(LSi == 0) then
   print *,'L = ',L,' atom = ',atom(L)
   print *,proximity(L,:)
   print *
   print *,'natom = ',natom-1,' atom = ',atom(natom-1)
   print *,proximity(natom-1,:)
   print *

return
end if
      r3v = rxyz(L,:) - rxyz(LSi,:)
      call pbc(r3v)
      r12 = r3v/len_3d(r3v)
      call zaxis2vect( r12, aa )
      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      zz = -bondl*cbeta
      xx = bondl*sbeta*cphi
      yy = bondl*sbeta*sphi
      r23 = matmul(aa,(/xx,yy,zz/))
      rxyz(LS,:) = rxyz(L,:) + r23
!
! attach the 3 Oxygen atoms
!
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      r23 = r23/len_3d(r23)
      call zaxis2vect(r23,aa)

      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      xx = rxy*sphi
      yy = -rxy*cphi
      tetra(1:3,1) = (/ xx, yy, zz /)
      tetra(1:3,2) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,3) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
!
      natom = natom + 1
      rxyz(natom,:) = rxyz(LS,:) + matmul(aa,tetra(1:3,1))
      natom = natom + 1
      rxyz(natom,:) = rxyz(LS,:) + matmul(aa,tetra(1:3,2))
      natom = natom + 1
      rxyz(natom,:) = rxyz(LS,:) + matmul(aa,tetra(1:3,3))
!
! Apply PBC to 4 added atoms
      call pbc(rxyz(natom-3:natom,:))
!
      call set_proximity(L,0,LS)
      atom(L) = iOxygen
      atom(LS) = iSilicon
      proximity(LS,1) = L
      proximity(LS,2) = LS + 1
      proximity(LS,3) = LS + 2
      proximity(LS,4) = LS + 3
      atom(LS+1:LS+3) = iOxygenH
      proximity(LS+1:LS+3,1) = LS
      atomL(natom-3:natom) = atomL(L)  ! cluster label for new atoms
!
      call NEW_NLIST(1,natom)
!
! Check for overlap
      overlap = .false.
      overlap_if: if (overlap_check) then
      outer_atom_loop: do ii = natom-3,natom
      ic = CELL(rxyz(ii,:))  ! link list
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if (atom(M) == iSilicon) GOTO 100
            if (M == ii) GOTO 100
            if (nearest_neighbor2(ii,m)) GOTO 100
            r3v(1:3) = rxyz(M,1:3) - rxyz(ii,1:3)
            call pbc(r3v)
            r2 = r3v(1)**2+r3v(2)**2+r3v(3)**2
            if (r2 <= (sigLJ_2(atom(m))+sigLJ_2(atom(ii)))**2) then
               overlap = .true.
               exit outer_atom_loop
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      end if overlap_if
!
      if (overlap) then
         call set_proximity(L,LS,0)
         atom(L) = iOxygenH
         proximity(LS:LS+3,:) = 0
         atom(LS:LS+3) = 0
         rxyz(LS:LS+3,:) = 0.0_wp
         natom = natom - 4
         call NEW_NLIST(1,natom)
         return
      else
         success = .true.
         return
      end if
   END SUBROUTINE attach_tetrahedra

END MODULE attachment_mod


!!>include 'water_attack_VL_mod.f90'

MODULE water_attack_mod
   USE precision_mod, only: wp
   implicit none
   real(wp):: fwater
CONTAINS

   SUBROUTINE water_attack(SiOH4_blocked,iat,success)
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE constants_mod, only: pi,angle_SiOH,bondl_oh
      USE verlet_list_mod
      USE nlist_mod
      USE rand_mod
      USE ran_point_sphere,only: ran3sph
      USE matvec3d_mod,only: len_3d,matvec_3d
      USE rotate_axis_mod
      USE HKNonLattice_mod
      logical,intent(in):: SiOH4_blocked
      integer,intent(in):: iat
      logical,intent(out):: success
      integer,parameter:: ntry=100
      integer:: copatom(natom_max),natom_old
      integer:: copproximity(natom_max,size(proximity,dim=2))
      real(wp):: coprxyz(natom_max,3)
      integer:: i,itmp,iSi1,iSi2
      real(wp):: r12(3) !!,sx,sy,sz,rad
      logical :: dangling_group
      success = .FALSE.
      dangling_group = .FALSE.
!
      ! The two attached Silicons
      iSi1 = proximity(iat,1)
      iSi2 = proximity(iat,2)
! if there is a dangling Si(OH)3 group the Silicon atom is iSi1
      if (is_dangling_group(iSi1,iat)) then
         dangling_group = .TRUE.
      end if
      if (is_dangling_group(iSi2,iat)) then
         ! if (dangling_group)STOP 'add_water: ERROR, both groups cannot be dangling'
         dangling_group = .TRUE.
         itmp = iSi2
         iSi2 = iSi1
         iSi1 = itmp
      end if
      if(dangling_group.and.SiOH4_blocked)then
         write(*,*) 'water attack & Si(OH)4 reattached'
         write(*,*) kmc_step,rxyz(iat,3)
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
! add the water
      atom(iat) = iOxygenH
      natom = natom + 1
!
!    O           OH     OH
!   / \   -->    |      |
! Si   Si        Si  +  Si
!
!    iat        natom       iat
!   /   \         |          |
! iSi1  iSi2     iSi1   +   iSi2
!
      atom(natom) = iOxygenH
      atomL(natom) = atomL(iSi1)   ! cluster label
      call set_proximity(natom,  0, iSi1)
      call set_proximity(iSi1, iat, natom)
      call set_proximity(iat, iSi1, 0)
! The oxygens overlap
      rxyz(natom,:) = rxyz(iat,:)
!
! move the Oxygens (iat and natom) to slightly better starting positions
      r12(1:3) = rxyz(iSi2,:) - rxyz(iSi1,:)
      call pbc(r12)
      r12 = r12/len_3d(r12)
      rxyz(iat,:) = rxyz(iat,:) + (0.1_wp/Angstrom)*r12
      rxyz(natom,:) = rxyz(natom,:) - (0.1_wp/Angstrom)*r12
!
!     call ran3sph(sx,sy,sz)
!     rad = sigLJ_2(atom(iat))
!     rxyz(natom,:) = rxyz(natom,:) + [sx,sy,sz]*rad
!
      call NEW_VLIST

1019  CONTINUE
      if (dangling_group) then
         call delete_group(iSi1)
         call NEW_VLIST
      end if
      write(*,*) 'add_water: success'
      write(79,*) kmc_step,rxyz(iat,3)
      success = .TRUE.
      RETURN

1020  CONTINUE ! the move failed so undo the changes
      rxyz(natom_old+1,:) = 0.0_wp
      proximity(natom_old+1,:) = 0
      atomL(natom_old+1) = 0
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
      call NEW_NLIST(1,natom)
   END SUBROUTINE water_attack

END MODULE water_attack_mod


!!>include 'condensation_mod.f90'

MODULE condensation_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE condensation(iO1,iO2)
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE verlet_list_mod
      USE nlist_mod
      integer,intent(inout):: iO1,iO2
      integer:: i,ii
      integer:: iSi2,iSi1
      real(wp):: r12(3)
!
!   Si-OH   +   Si-OH  -->  Si-O-Si  +  H2O
!
!     iO1   +  iO2           iO1           iO2
!     /        /     -->    /    \    +   /  \
!   iSi1     iSi2         iSi1  iSi2
!

!     Pick iSi1, the Silicon attached to iO1
      iSi1 = 0
      do i = 1,ncmax(atom(iO1))
         ii = proximity(iO1,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi1 = ii
         end if
      end do
      iSi2 = 0
      do i = 1,ncmax(atom(iO2))
         ii = proximity(iO2,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi2 = ii
         end if
      end do

! new position for iO1
      r12 = rxyz(iO2,:) - rxyz(iO1,:)
      call pbc(r12)
      rxyz(iO1,:) = rxyz(iO1,:) + 0.5_wp*r12
      call pbc(rxyz(iO1,:))
!
      atom(iO1) = iOxygen
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO1, 0, iSi2)
      call delete_atom(iO2)
!
      call new_vlist
!
      write(*,*) 'H2O removal ',rxyz(iO1,3)
      write(80,*) kmc_step,rxyz(iO1,3)
   END SUBROUTINE condensation

END MODULE condensation_mod


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
      USE global_vars_mod, only: natom_max
      USE nlist_mod
      USE atom_types_mod
      USE constants_mod, only: sigma_2
      USE ran_point_sphere,only: ran3sph
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: rprobe
      integer:: ia,it,ic,j,jj,nc,neigh,ncell(125)
      real(wp):: sx,sy,sz,r0(3)
      real(wp):: dr(3),r(3),rad
      if (.not.allocated(surface)) allocate(surface(natom_max))
      do ia = ibegin,iend
         surface(ia) = .false.
if (atom(ia) == iSilicon) cycle
         rad = rprobe + sigLJ_2(atom(ia))
         r0(1:3) = rxyz(ia,1:3)
         trial_loop: do it = 1,ntrial
            call ran3sph(sx,sy,sz)
            r(1:3) = r0(1:3) + (/sx,sy,sz/)*rad
            ic = CELL(r)
            CALL NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj = 1,neigh
               nc = ncell(jj)
               if (nc == 0)cycle cell_loop
               j = HOC(nc)
               do while (j /= 0)
if (atom(j) == iSilicon) GOTO 100
                  if (j == ia) GOTO 100
                  dr(1:3) = r(1:3) - rxyz(j,1:3)
                  call pbc(dr)
                  if (dot_product(dr,dr) < (rprobe + sigLJ_2(atom(j)))**2) then
                     cycle trial_loop
                  end if
100               j = LL(j)
               end do
            end do cell_loop
            surface(ia) = .true.
            exit trial_loop
         end do trial_loop
      end do
      fsurfa = real(count(surface(ibegin:iend)),wp)/(iend-ibegin+1)
   END SUBROUTINE

   SUBROUTINE write_surface_atoms_zbin(iu)
      USE deposition_mod
      USE nlist_mod
      integer,intent(in):: iu
      integer:: iz,ic,nsur,nat,nc,L,neigh,ncell(125)
      real(wp):: fsur
      do iz = 1,nbin
      nsur = 0
      nat = 0
      CALL Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
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
      write(iu,'(5f12.6)') (iz-0.5_wp)*delz,fsur,fsur/voidage(iz)
      end do
   END SUBROUTINE

END MODULE SURFACE_ATOMS_MOD

!!>include 'switch.f90'
MODULE bond_switch_O_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE bond_switch_O(success)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      USE global_vars_mod
      USE nlist_mod
      logical,intent(out):: success
      integer:: iOb,iO1,iO2
      integer:: iSi1,iSi2,iSi3,iSi4
!      logical:: consistent

! TGV method
! Bond switch for SiO2 proposed by:
! Yuhai Tu, G. Grinstein, and D. Vanderbilt, Phys. Rev. Lett. 81,4899 (1998).
!
!      O        O2              O  O2
!      |        |               | /
!  O--Si1--Ob--Si2--O  -->  O--Si1--Ob--Si2--O
!      |        |                      / |
!      O1       O                    O1  O
!

! randomly select a bridging oxygen iOb
      do
         iOb = int(rand()*natom) + 1
         if (atom(iOb) == iOxygen) exit
      end do
      iSi1 = proximity(iOb,1)
      iSi2 = proximity(iOb,2)
      if (iSi1==0) return ! stop 'oxygen malfunction 1'
      if (iSi2==0) return ! stop 'oxygen malfunction 2'

! randomly select an oxygen iO1 attached to iSi1 (not iOb)
      do
         iO1 = proximity(iSi1,int(rand()*ncmax(atom(iSi1)))+1)
         if (iO1 /= iOb) exit
      end do
! repeat for iSi2
      do
         iO2 = proximity(iSi2,int(rand()*ncmax(atom(iSi2)))+1)
         if (iO2 /= iOb) exit
      end do
if (iO1 == 0) return
if (iO2 == 0) return
! no 2-Si rings are allowed to form
      iSi3 = proximity(iO2,1)
      if (iSi3 == iSi2) iSi3 = proximity(iO2,2)
      if (nearest_neighbor2(iSi3,iSi1)) then
         success = .FALSE.
         RETURN
      end if
      iSi4 = proximity(iO1,1)
      if (iSi4 == iSi1) iSi4 = proximity(iO1,2)
      if (nearest_neighbor2(iSi4,iSi2)) then
         success = .FALSE.
         RETURN
      end if

      call set_proximity(iSi1,iO1,iO2)
      call set_proximity(iO1,iSi1,iSi2)
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO2,iSi2,iSi1)
      success = .TRUE.

   END SUBROUTINE bond_switch_O

END MODULE bond_switch_O_mod

!!>include 'relax_bonds_VL_NL_mod.f90'

MODULE relax_bonds_mod
   USE precision_mod, only: wp
   USE bond_switch_O_mod
   implicit none
CONTAINS

   SUBROUTINE relax_bonds(nouter)
      USE global_vars_mod
      USE constants_mod, only: bondl_oh
      USE coordinates_mod
      USE connectivity_mod
      USE Relax_mod
      USE Keating
      USE repul_energy_mod
      USE rand_mod
      USE HKNonLattice_mod
      USE verlet_list_mod
      USE nlist_mod
      USE sort_mod
      integer,intent(in):: nouter
      integer:: copproximity(natom_max,size(proximity,dim=2)),copatom(natom_max)
      real(wp):: coprxyz(natom_max,3)
      real(wp):: EEE,energy1,energy2
      integer:: iouter,i,natom_old
      real(wp):: del_energy
      logical :: success
!
      outer_loop: do iouter = 1,nouter

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

      call bond_switch_O(success)
      if (.NOT.success)then
         write(*,*) 'bond switch ',' failed: connectivity'
         cycle
      end if
!
      call RELAX(NRELAX, 1, natom, kboltzT)
!
! energy after the move
      energy2 = energy4(1,natom) + repul_energy(1,natom)
      del_energy = energy2 - energy1
      EEE = del_energy/kboltzT
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
      call NEW_NLIST(1,natom)
      write(*,*) 'bond switch ',' failed: energy'
      cycle outer_loop

1019  CONTINUE
      write(*,*) 'bond switch ','success' !,rxyz(iat,3)

      end do outer_loop
      RETURN
   END SUBROUTINE

END MODULE relax_bonds_mod
!!>include 'find_atoms_mod.f90'

MODULE find_atoms_mod
   implicit none
   public:: find_atoms
CONTAINS
   SUBROUTINE find_atoms(iat,atomtype,rad,nat,atomlst)
      use nlist_mod
      use atom_types_mod
      use coordinates_mod
      integer,intent(in):: iat,atomtype
      real(wp),intent(in):: rad
      integer,intent(out):: nat,atomlst(:)
      integer:: ic,nc,j,jj,nl,neigh,ncell(125)
      real(wp):: dr(3),r2
! Look for atoms of type atomtype within a distance rad of ia
      r2 = rad*rad
      nat = 0
      atomlst = 0
      nl = nlayers_ll(rad)
      ic = CELL(rxyz(iat,:))
      call NEIGCELL(ic,nl,neigh,ncell)
      cell_loop: do jj = 1,neigh
      nc = ncell(jj)
      if (nc == 0)cycle cell_loop
      j = HOC(nc)
      do while (j /= 0)
         if (j == iat) GOTO 100
         if (atom(j) /= atomtype) GOTO 100
         dr(1:3) = rxyz(j,1:3) - rxyz(iat,1:3)
         call pbc(dr(:))
         if ((dr(1)**2+dr(2)**2+dr(3)**2) <= r2) then
            nat = nat + 1
            atomlst(nat) = j
         end if
100      j = LL(j)
      end do
      end do cell_loop
   END SUBROUTINE

END MODULE

!!>include 'list_atoms_in_cells_mod.f90'

MODULE list_atoms_in_cells_mod
   USE precision_mod
   implicit none
   public:: list_atoms_in_cells
   integer,parameter:: L_OH = 1
   integer,parameter:: L_OSiOH3 = 2
   integer,parameter:: L_OSiO = 3
!
   integer,parameter:: iOSiOH3 = 6
   integer,parameter:: iOSiO = 7
   integer,parameter:: iOHB = 8
   integer,parameter,private:: nbinmax = 1000
   real(wp):: volbin
   real(wp):: ndenrxn(-nbinmax:nbinmax,0:8)

CONTAINS

   SUBROUTINE list_atoms_in_cells(numcells,naborcells,alst)
      USE nlist_mod
      USE list_mod
      USE atom_types_mod
      USE connectivity_mod
      implicit none
      integer,intent(in):: numcells
      integer,intent(in):: naborcells(:)
      type(list),intent(out):: alst(:)
      integer:: i,ia,ic,nc,L
      alst(:)%n = 0
      forall(i=1:size(alst)) alst(i)%ind = 0
      cell_loop: do ic = 1,numcells
         nc = naborcells(ic)
         if (nc == 0)cycle cell_loop
         ia = HOC(nc)
         atom_loop: do while (ia /= 0)
            select case(atom(ia))
            case(iOxygenH)
               L = L_OH
            case(iOxygen)
               if ( is_dangling_group(proximity(ia,1),ia).or. &
                    is_dangling_group(proximity(ia,2),ia) )then
                  L = L_OSiOH3
               else
                  L = L_OSiO
               end if
            case default
               goto 100
            end select
            alst(L)%n = alst(L)%n + 1
            alst(L)%ind(alst(L)%n) = ia
100         ia = LL(ia)
         end do atom_loop
      end do cell_loop
   END SUBROUTINE


   SUBROUTINE den_profile_rxn(ibegin,iend,dl)
      USE coordinates_mod
      USE atom_types_mod
      USE connectivity_mod
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,dli
      integer:: j,ii,jj,iSi1,iSi2,nbin
      top_atom = maxval(rxyz(ibegin:iend,3))
      nbin = int(top_atom/dl) + 1
      dli = 1.0_wp/dl
      ! volbin = boxl(1)*boxl(2)*dl
      ndenrxn(:,:) = 0.0_wp
      do j = ibegin,iend
         ii = int(rxyz(j,3)*dli + 1.0_wp)
         if (ii > nbin) then
            if ( rxyz(j,3) == top_atom) then
               ii = nbin
            else
               stop 'ii > nbin'
            end if
         end if
         jj = atom(j)
         ndenrxn(ii,jj) = ndenrxn(ii,jj) + 1.0_wp
         if(jj == iOxygen)then
            iSi1 = proximity(j,1)
            iSi2 = proximity(j,2)
            if(iSi1 == 0) iSi1 = iSi2
            if (is_dangling_group(iSi1,j).or.is_dangling_group(iSi2,j)) then
               ndenrxn(ii,iOSiOH3) = ndenrxn(ii,iOSiOH3) + 1.0_wp
            else
               ndenrxn(ii,iOSiO) = ndenrxn(ii,iOSiO) + 1.0_wp
            end if
         end if
      end do
      ! ndenrxn(1:nbin,:) = ndenrxn(1:nbin,:)/volbin
   END SUBROUTINE

END MODULE



PROGRAM porous_silica
      USE precision_mod, only: wp
      USE global_vars_mod
      USE constants_mod
      USE atom_types_mod
      USE list_mod
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      USE frames_mod
      USE Relax_mod
      USE nlist_mod
      USE verlet_list_mod
      USE HKNonLattice_mod
      USE density_profile_mod
      USE deposition_mod
      USE rdf_mod
      USE Keating_parameters, only: ASiO
      USE water_attack_mod
      USE files_mod
      USE Relax_mod
      USE select_atom_mod
      USE attachment_mod
      USE probe_mol_mod
      USE voidage_mod
      USE surface_atoms_mod
      USE keating,only:energy4
      USE repul_energy_mod
      USE seaton_mod
      USE charges_mod
      USE lj_el_mod
      USE Henrys_law_calc_mod
      USE condensation_mod
      USE relax_bonds_mod
      USE find_atoms_mod
      USE list_atoms_in_cells_mod
      implicit none
      real(wp),parameter:: ulength = 1.0e-10_wp*angstrom
      real(wp),parameter:: uenergy = 1.60217733e-19_wp
      real(wp),parameter:: upressure = uenergy/(ulength**3)
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: aSiOSi = pi*(140.0_wp/180.0_wp)
      real(wp):: top_atom,void_crit,rverlet,dlayer
      integer:: i,itmp,j,ipconfig,ntotal,imve,narg
      real:: timep(0:10)
      character:: date*8,ctime*10
      character(len=32):: ctmp,c6*6,c5*5
      character(len=132):: carg
      logical:: success,ok,booltmp
!
! Kinetic Monte Carlo
integer,parameter:: nrxn = 4
integer,parameter:: nbinmax = 10000
real(wp):: rate(nrxn,nbinmax),rate_sum,psum
integer:: irxn(nrxn,nbinmax,0:3)
type(list):: alst(3),oh_lst,ohb_lst
integer:: ia_OH
integer:: ia_OSiOH3
integer:: ia_OSiO
integer:: ib,ir,iOH,len,status
real(wp):: time,delta_time,xi1
real(wp),parameter:: rohbl = 4.0_wp/Angstrom
real(wp),parameter:: volohbl = (4.0_wp/3.0_wp)*pi*(rohbl**3)
real(wp),parameter:: torr2Pa = 101325.0_wp/760.0_wp
integer:: ntmp,neigh,ncell(125)
real(wp):: P_SiOH4,P_H2O,exp_Ea_kt,Temp_K,zu,zl,vco2,vo2,vn2
logical:: hydrothermal = .false.
!      OPEN(unit = 10,FILE = 'kmc.dat',STATUS = 'OLD')
!-----------------------------------------------------------
!     Read in simulation data
!     kmc_step = number of kmc steps so far
!     ntotal = total number of kmc steps
!     boxl = box side length of cell
!     ipconfig = frequency of xyz,RDF dump
!     irand_seed = random number seed
!-----------------------------------------------------------
!      READ(10,*) kmc_step
!      READ(10,*) ntotal
!      READ(10,*) Temp_K
!      READ(10,*) P_SiOH4
!      READ(10,*) fwater
!      READ(10,*) boxl
!      READ(10,*) void_crit
!      READ(10,*) ipconfig
!      READ(10,*) irand_seed
!      READ(10,*) NRELAX
!      READ(10,*) rverlet
!      READ(10,*) ntmp
!      READ(10,*) e_activ,del_rxn

      narg = command_argument_count()
      call get_command_argument (0, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', carg(1:len)
      if(narg < 5)then
         write(*,*)'usage :'
         write(*,*) carg(1:len),' kmc_step  ntotal  T[K]   P_H2O/P_Si(OH)4   random_seed'
         stop
      end if
!
      call get_command_argument (1, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 1
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) kmc_step
      print *,'kmc_step = ',kmc_step
!
      call get_command_argument (2, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 1
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) ntotal
      print *,'ntotal = ',ntotal

      call get_command_argument (3, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 3
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) Temp_K

      call get_command_argument (4, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 4
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) fwater

      call get_command_argument (5, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 5
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) irand_seed
      print *,'irand_seed = ',irand_seed


      P_SiOH4 = 1333.22375_wp
      void_crit=0.035
      ipconfig=50
      NRELAX=10
      rverlet=0.5
      ntmp=1
      e_activ=0.75;del_rxn=0.14
!      close(UNIT = 10,STATUS = 'KEEP')
!
      kboltzT = K_ev*Temp_K
      exp_Ea_kt = exp(-del_rxn/kboltzT)
      P_SiOH4 = P_SiOH4/upressure
      P_H2O = fwater*P_SiOH4
if(hydrothermal) then
   P_H2O = P_SiOH4
   P_SiOH4 = 0.0_wp
end if
!
      print '(a,g18.9,a)','ulength   = ',ulength,' m'
      print '(a,g18.9,a)','uenergy   = ',uenergy,' J'
      print '(a,g20.10,a)','upressure = ',upressure,' Pa'
      print '(a,g18.9,a,g18.9,a)','P_SiOH4 = ',P_SiOH4,' = ',P_SiOH4*upressure,' Pa'
      print '(a,g18.9,a)',        '        = ',P_SiOH4*upressure/torr2Pa,' torr'
      print '(a,g18.9,a)','kT   = ',kboltzT,' eV'
      print *,'T = ',Temp_K,' K'
      print '(a,g18.9)','exp(-EA/kT) = ',exp_Ea_kt
!
      natom_max = 10000 + 7*ntotal
      allocate(rxyz(1:natom_max,3))
      allocate(atom(0:natom_max)); atom(0) = -1
      allocate(charge(1:natom_max))
      allocate( proximity(natom_max,4) )
      proximity = 0
      call Init_Verlet_List(rverlet,0.26_wp)
      call Init_HKNonLattice(natom_max)
      call LJ_init
!
      OPEN(unit= 15,file= 'psil_config.pdb')
      OPEN(unit= 20,file= 'psil_topology.out')
      OPEN(unit= 14,file= 'psil_RDF.out',POSITION='APPEND')
      OPEN(unit= 21,file= 'psil_rand.out')
      OPEN(unit= 22,file= 'psil_dens.out')
      OPEN(unit= 23,file= 'psil_void_SiO4.out')
      OPEN(unit= 26,file= 'psil_void_H2O.out')
      OPEN(unit= 24,file= 'psil_topatom.out')
      OPEN(unit= 25,file= 'psil_time.out')
      OPEN(unit= 29,file= 'psil_kmc_history.out')
      OPEN(unit= 77,file= 'psil_SiO4_attempt.out')
      OPEN(unit= 78,file= 'psil_SiO4_add.out')
      OPEN(unit= 79,file= 'psil_H2O_add.out')
      OPEN(unit= 80,file= 'psil_H2O_remove.out')
      OPEN(unit= 81,file= 'psil_bond_switch.out')


!      call read_xyz(15,kmc_step)
!      call read_proximity(20,kmc_step)
      read (15,*) c6,natom  ;print *,'natom = ',natom
      read (15,*) c6,boxl   ;print '(a,6f14.8)',' boxl = ',boxl
      do i = 1,natom
         read (15,'(a6)',advance = 'no') c6
         read (15,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
!print '(i6,3f12.6,i6)',i,rxyz(i,:),atom(i)
      end do
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
      do i = 1,natom
         read(15,'(a32)') ctmp
         do j = 1,ncmax(atom(i))
            c5 = ctmp(6 + 5*j + 1:6 + 5*j + 5)
            read( unit=c5,fmt=* ) proximity(i,j)
! print *,i,' -->  ',proximity(i,:)
         end do
      end do
!
      boxl = boxl/Angstrom
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!     boxl2n = -boxl2
!     boxl2i = 1.0_wp/boxl2

!
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
      write(*,*) 'n_cluster = ',n_cluster
      write(*,*) 'atomL = '; write(*,'(10i6)') atomL(1:natom)
      top_atom = maxval(rxyz(1:natom,3))
!
      call init_probe_mols()
      call INIT_NLIST(boxl(1),boxl(2),500.0_wp/angstrom,max(5.0_wp/angstrom,2.0_wp*maxval(probe_tip3p%rad)),natom_max)
      call NEW_NLIST(1,natom)
      call new_vlist
      call cpu_time(timep(0))
!
!----------------------------------------------------------------------
      time = 0.0_wp
      main_loop: do while (kmc_step < ntotal)
      call cpu_time(timep(1))
      kmc_step = kmc_step+1
!
      call assign_charge(1,natom)

      dlayer = (5.0_wp/angstrom)
      volbin = boxl(1)*boxl(2)*dlayer
!      call henry_profile_z(probe_SiO4,dlayer,1)
!      call henry_profile_z(probe_tip3p,dlayer,2)
      call henry_profile(probe_SiO4,dlayer,1)
      call henry_profile(probe_tip3p,dlayer,2)
      call den_profile_rxn(1,natom,dlayer)
      call prob_dist(henryk(1:nbin,1),inthenryk(1:nbin,1),void_crit,itop_nonp_SiO4)
      call prob_dist(henryk(1:nbin,2),inthenryk(1:nbin,2),void_crit,itop_nonp_H2O)
!do i = 1,nbin
!   write(99,'(i5,2(1x,f0.6))')i,henryk(i,1),henryk(i,2)
!end do
!write(99,*)
!
!                         k_1
! (1)  |-OH  +  Si(OH)4   --->   |-Si(OH)3  +  H2O
!
!                         k_1R
! (2)  |-Si(OH)3  +  H2O  --->   |-OH  +  Si(OH)4
!
!                         k_2
! (3)     |-OH  +  HO-|   --->   |-O-|  +  H2O
!
!                         k_2R
! (4)     |-O-|  +  H2O   --->   |-OH  +  HO-|
!
      rate_sum = 0.0_wp
      rate(1:nrxn,1:nbin) = 0.0_wp
      irxn(1:nrxn,1:nbin,:) = 0
      do ib = 1,nbin
         call Z_NEIGH_CELL(ib,neigh,ncell)
         call list_atoms_in_cells(neigh,ncell,alst)
         call rand_from_list(alst(L_OH),ia_OH)
         call rand_from_list(alst(L_OSiOH3),ia_OSiOH3)
         call rand_from_list(alst(L_OSiO),ia_OSiO)
         if(ia_OSiO /= 0)then
            rate(4,ib) = exp_Ea_kt*henryk(ib,2)*P_H2O*ndenrxn(ib,iOSiO)
            if(ib < itop_nonp_H2O) rate(4,ib)=0.0_wp
            rate_sum = rate_sum + rate(4,ib)
            irxn(4,ib,0) = 1
            irxn(4,ib,1) = ia_OSiO
         end if
         if(ia_OSiOH3 /= 0)then
            rate(2,ib) = exp_Ea_kt*henryk(ib,2)*P_H2O*ndenrxn(ib,iOSiOH3)
            if(ib < itop_nonp_H2O) rate(2,ib) = 0.0_wp
            rate_sum = rate_sum + rate(2,ib)
            irxn(2,ib,0) = 1
            irxn(2,ib,1) = ia_OSiOH3
         end if
         if(ia_OH /= 0)then
            rate(1,ib) = henryk(ib,1)*ndenrxn(ib,iOxygenH)*P_SiOH4
            if(ib < itop_nonp_SiO4) rate(1,ib) = 0.0_wp
            rate_sum = rate_sum + rate(1,ib)
            irxn(1,ib,0) = 1
            irxn(1,ib,1) = ia_OH
            call find_atoms(ia_OH, iOxygenH, rohbl, oh_lst%n, oh_lst%ind)
            ohb_lst%n = 0
            ohb_lst%ind = 0
            do i = 1,oh_lst%n
               iOH = oh_lst%ind(i)
               if(OH_groups_OK(ia_OH,iOH)) call add_to_list(ohb_lst,iOH)
            end do
            call rand_from_list(ohb_lst,iOH)
            if(iOH /= 0)then
               rate(3,ib) = kboltzT*(ndenrxn(ib,iOxygenH)/volbin)*ohb_lst%n*(volbin/volohbl)
               if(ib <  itop_nonp_H2O) rate(3,ib) = 0.0_wp
               rate_sum = rate_sum + rate(3,ib)
               irxn(3,ib,0) = 2
               irxn(3,ib,1) = ia_OH
               irxn(3,ib,2) = iOH
            end if
         end if
      end do

      xi1 = rand()
      psum = 0.0_wp
      bin_loop: do ib = 1,nbin
      do ir = 1,nrxn
         psum = psum + rate(ir,ib)/rate_sum
         if (psum >= xi1) then
            exit bin_loop
         end if
      end do
      end do bin_loop
      write(*,'(g18.8,i4,2i9)') rate_sum,ir,ib,nbin
      write(29,'(g18.8,i4,2i9)') rate_sum,ir,ib,nbin
!
!     carry out reaction #ir in bin #ib
print *,'ir = ',ir

do
   ia_OH = int(rand()*natom) + 1
   if (atom(ia_OH) == iOxygenH) exit
end do
call find_atoms(ia_OH, iOxygenH, rohbl, oh_lst%n, oh_lst%ind)
call rand_from_list(ohb_lst,iOH)
if (ir > 4)then
if(rand()<0.9)then
    ir = 1
else
    ir = 3
  if(iOH==0)ir=1
end if
end if
print *,'ir = ',ir
      select case(ir)
      case(1) ! reaction 1
!         ia_OH = irxn(1,ib,1)
         call attach_tetrahedra(ia_OH,aSiOSi,aSiO,overlap_check=.false.,success=success)
print *,'##################### 1 ###############'
         if(success)then
            write(77,*) kmc_step,ib,rxyz(ia_OH,3)
         end if
      case(2) ! reaction 2
         ia_OSiOH3 = irxn(2,ib,1)
         booltmp = (ib < itop_nonp_SiO4)
         call water_attack(booltmp,ia_OSiOH3,success)
print *,'##################### 2 ###############'
      case(3)
!         ia_OH = irxn(3,ib,1)
!         iOH  =  irxn(3,ib,2)
         call condensation(ia_OH,iOH)
print *,'##################### 3 ###############'
      case(4)
         ia_OSiO = irxn(4,ib,1)
         booltmp = (ib < itop_nonp_SiO4)
         call water_attack(booltmp,ia_OSiO,success)
print *,'##################### 4 ###############'
      end select
!
!     advance the kmc clock
      delta_time = -log(rand())/rate_sum
      time = time + delta_time
!
      call cpu_time(timep(2))
      call new_vlist
!
      call RELAX(NRELAX, 1, natom, kboltzT)

      call cpu_time(timep(3))

!     call relax_bonds( count(atom(1:natom) == iOxygen) )
      call relax_bonds( min(10,count(atom(1:natom) == iOxygen)) )


      call cpu_time(timep(4))
      call date_and_time(date,ctime)

!     write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
!     write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
      write(* ,'(i7,a10,2x,a10,6(3x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
      write(25,'(i7,a10,2x,a10,6(2x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
      top_atom = maxval(rxyz(1:natom,3))
      write(24,'(f12.6,2i7,f12.6)') time,kmc_step,natom,top_atom
!
      if (make_movie) then
         call new_frame_file(imve,'frame',kmc_step)
         call write_frame(imve,1,natom)
      end if
!
      if ( mod(kmc_step,ipconfig) == 0 .or. (kmc_step == ntotal) ) then
!--------Periodically check the connectivity and write out xyz coordinates,
!--------RDF, Connectivity, etc.
!         call check_proximity(ok,i)
!         if (.NOT.ok) then
!            write(*,*) 'proximity array is NOT consistent'
!            write(*,*) 'check_proximity ',ok,i
!            stop
!         end if
!         call check_proximity2(1,natom,ok,i)
!         if (.NOT.ok) then
!           write(*,*) 'proximity array is NOT consistent'
!           write(*,*) 'check_proximity2 ',ok,i
!           stop
!         end if
         call rdf_calc
         call rdf_print(14,kmc_step)
         call write_xyz(15)
         rewind(15)
         call write_proximity(20,kmc_step)
         rewind(20)
         call write_rand_state(21)
         rewind(21)
         call den_profile(1,natom,(5.0_wp/angstrom))
         call write_den_profile_z(22,kmc_step)
!         call voidage_profile(probe_SiO4,dlayer)
!         call write_voidage(23)
!         call SURFACE_ATOMS(1,natom,(1.5_wp/angstrom))
!         write(*,*) 'fraction of surface atoms = ',fsurfa
!         !call write_surface_atoms_zbin(191)
!         call voidage_profile(probe_tip3p,dlayer)
!         call write_voidage(26)
      end if

      end do main_loop
!----------------------------------------------------------------------




!
      write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
      write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
      write(*,*) 'N_OH= ',count(atom(1:natom) == iOxygenH)
      write(*,*) 'energy  ', energy4(1,natom), repul_energy(1,natom)

!
      close(14)
      OPEN(unit= 14,file= 'analysis_RDF.out')
      call rdf_calc
      call rdf_print(14,kmc_step)
      close(14)

      close(22)
      OPEN(unit= 22,file= 'analysis_dens.out')
      call den_profile(1,natom,(10.0_wp/angstrom))
      call write_den_profile_z(22,kmc_step)
      close(22)

!      close(23)
!      OPEN(unit= 23,file= 'analysis_void_SiO4.out')
!      call voidage_profile(probe_SiO4,dlayer)
!      call write_voidage(23)
!      close(23)
!
!
!      close(24)
!      OPEN(unit= 24,file= 'analysis_void_H2O.out')
!      call voidage_profile(probe_H2O,dlayer)
!      call write_voidage(24)
!      close(24)
!
!      close(25)
!      OPEN(unit= 25,file= 'analysis_void_CO2.out')
!      call voidage_profile(probe_CO2,dlayer)
!      call write_voidage(25)
!      close(25)
!
!      close(26)
!      OPEN(unit= 26,file= 'analysis_void_N2.out')
!      call voidage_profile(probe_N2,dlayer)
!      call write_voidage(26)
!      close(26)
!
!      close(27)
!      OPEN(unit= 27,file= 'analysis_void_O2.out')
!      call voidage_profile(probe_O2,dlayer)
!      call write_voidage(27)
!      close(27)

      zl =  5.0_wp/angstrom
      zu = 25.0_wp/angstrom
      vo2 = voidage_calc(50000,zl,zu,probe_O2)
      vn2 = voidage_calc(50000,zl,zu,probe_N2)
      vco2 = voidage_calc(50000,zl,zu,probe_CO2)
      open(unit=40,file='analysis_void_5-25.out')
      write(*,*) zl,zu
      write(*,*) 'voidage O2 ',vO2
      write(*,*) 'voidage N2 ',vN2
      write(*,*) 'voidage CO2 ',vCO2
      write(40,*) zl,zu
      write(40,*) 'voidage O2 ',vO2
      write(40,*) 'voidage N2 ',vN2
      write(40,*) 'voidage CO2 ',vCO2
      close(40)
!
      zl =  5.0_wp/angstrom
      zu = 15.0_wp/angstrom
      vo2 = voidage_calc(50000,zl,zu,probe_O2)
      vn2 = voidage_calc(50000,zl,zu,probe_N2)
      vco2 = voidage_calc(50000,zl,zu,probe_CO2)
      open(unit=41,file='analysis_void_5-25.out')
      write(*,*) zl,zu
      write(*,*) 'voidage O2 ',vO2
      write(*,*) 'voidage N2 ',vN2
      write(*,*) 'voidage CO2 ',vCO2
      write(41,*) zl,zu
      write(41,*) 'voidage O2 ',vO2
      write(41,*) 'voidage N2 ',vN2
      write(41,*) 'voidage CO2 ',vCO2
      close(41)

!

!      close(28)
!      OPEN(unit= 28,file= 'analysis_henry_SiO4.out')
!      call henryslaw_profile(probe_SiO4,dlayer)
!      call prob_dist(void_crit)
!      itop_nonp_SiO4 = itop_nonp
!      call write_voidage(28)
!      call write_voidage_dist(28)
!      close(28)
!
!      close(29)
!      OPEN(unit= 29,file= 'analysis_henry_SiOH4.out')
!      call henryslaw_profile(probe_SiOH4,dlayer)
!      call prob_dist(void_crit)
!      itop_nonp_SiO4 = itop_nonp
!      call write_voidage(29)
!      call write_voidage_dist(29)
!      close(29)
!
      close(30)
      OPEN(unit= 30,file= 'analysis_surface.out')
      call SURFACE_ATOMS(1,natom,(1.5_wp/angstrom))
      write(*,*) 'fraction of surface atoms = ',fsurfa
      call write_surface_atoms_zbin(30)
      close(30)
!
      STOP
END PROGRAM porous_silica

