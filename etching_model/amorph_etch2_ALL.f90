!!>INCLUDE 'precision_mod.f90'

MODULE precision_mod
   implicit none
!  integer, parameter:: sp = kind(1.0)
!  integer, parameter:: dp = kind(1.0d0)
   integer, parameter:: sp = selected_real_kind(6,30)
   integer, parameter:: dp = selected_real_kind(15,300)
   integer, parameter:: qp_preferred = selected_real_kind(30,1000)
   integer, parameter:: qp = (1 + sign(1,qp_preferred))/2*qp_preferred &
                           + (1 - sign(1,qp_preferred))/2*dp
!
   integer,parameter,public:: wp = dp
   integer,parameter,public:: i2b = selected_int_kind(4)
   integer,parameter,public:: i4b = selected_int_kind(9)
   integer,parameter,public:: i8b = selected_int_kind(18)
END MODULE precision_mod

!!>INCLUDE 'math_const_mod.f90'

MODULE math_const_mod
!!
!!***********************************************************************
!!
!!     MODULE : MATH_CONST_MOD
!!     ACTION : Define some mathematical constants
!!     SOURCE : Calculated using Mathematica(R) 4.1
!!              http://www.wolfram.com
!!
!!       $Author: TmcD $
!!       $Date: 2001/09/17 09:41:40 $
!!       $Revision: 1.0 $
!!
!!***********************************************************************
!!
      USE precision_mod, only: wp
      implicit none
      private:: wp
      public
!
      real(wp), parameter:: m_e        = 2.71828182845904523536028747135266249775725_wp  ! e
      real(wp), parameter:: m_log2e    = 1.44269504088896340735992468100189213742665_wp  ! log_2 e
      real(wp), parameter:: m_log10e   = 0.434294481903251827651128918916605082294397_wp ! log_10 e
      real(wp), parameter:: m_ln2      = 0.693147180559945309417232121458176568075500_wp ! log_e 2
      real(wp), parameter:: m_ln10     = 2.30258509299404568401799145468436420760110_wp  ! log_e 10
      real(wp), parameter:: m_pi       = 3.14159265358979323846264338327950288419717_wp  ! pi
      real(wp), parameter:: m_pi_2     = 1.57079632679489661923132169163975144209858_wp  ! pi/2
      real(wp), parameter:: m_pi_4     = 0.785398163397448309615660845819875721049292_wp ! pi/4
      real(wp), parameter:: m_1_pi     = 0.318309886183790671537767526745028724068919_wp ! 1/pi
      real(wp), parameter:: m_2_pi     = 0.636619772367581343075535053490057448137839_wp ! 2/pi
      real(wp), parameter:: m_sqrtpi   = 1.77245385090551602729816748334114518279755_wp  ! sqrt(pi)
      real(wp), parameter:: m_1_sqrtpi = 0.564189583547756286948079451560772585844051_wp ! 1/sqrt(pi)
      real(wp), parameter:: m_2_sqrtpi = 1.12837916709551257389615890312154517168810_wp  ! 2/sqrt(pi)
      real(wp), parameter:: m_sqrt2    = 1.41421356237309504880168872420969807856967_wp  ! sqrt(2)
      real(wp), parameter:: m_sqrt3    = 1.73205080756887729352744634150587236694281_wp  ! sqrt(3)
      real(wp), parameter:: m_sqrt1_2  = 0.707106781186547524400844362104849039284836_wp ! 1/sqrt(2)
!
      real(wp), parameter:: m_euler    = 0.577215664901532860606512090082402431042159_wp   ! euler gamma
      real(wp), parameter:: m_2pi      = 6.28318530717958647692528676655900576839434_wp    ! 2*pi
      real(wp), parameter:: m_lnpi     = 1.14472988584940017414342735135305871164729_wp    ! log_e pi
      real(wp), parameter:: m_deg2rad  = 0.0174532925199432957692369076848861271344287_wp  ! deg2rad
      real(wp), parameter:: m_rad2deg  = 57.2957795130823208767981548141051703324055_wp    ! rad2deg
!
END MODULE math_const_mod

!!>INCLUDE 'phys_const_mod.f90'

MODULE phys_const_mod
!!
!!***********************************************************************
!!
!!     MODULE : PHYS_CONST_MOD
!!     ACTION : Define some physical constants
!!     SOURCE : Fundamental Physical Constants from NIST (1998 values)
!!              http://physics.nist.gov/cuu/Constants/
!!
!!       $Author: TmcD $
!!       $Date: 2001/09/17 09:51:35 $
!!       $Revision: 1.0 $
!!
!!***********************************************************************
!!
      USE precision_mod, only: wp
      USE math_const_mod, only: m_pi,m_2pi
      implicit none
      private:: wp,m_pi,m_2pi
      public
!
      real(wp),parameter:: amu       =  1.66053873e-27_wp      ! u        kg
      real(wp),parameter:: c_light   =  299792458.0_wp         ! c        m*s(-1)
      real(wp),parameter:: c2_light  =  c_light*c_light        ! c^2      m(2)*s(-2)
      real(wp),parameter:: navo     =  6.02214199e+23_wp       ! n_a      mol(-1)
      real(wp),parameter:: eplus    =  1.602176462e-19_wp      ! e        C
      real(wp),parameter:: eplus2   =  eplus*eplus             ! e^2      C(2)
      real(wp),parameter:: mu0      =  m_pi*4.0e-7_wp          ! mu_0     n*a(-2)
      real(wp),parameter:: eps0 =  1.0_wp/(c2_light*mu0)       ! eps_0    f*m(-1)
      real(wp),parameter:: elm_coupling = eplus2/(4.0_wp*m_pi*eps0)
      real(wp),parameter:: melectron  =  9.10938188e-31_wp    ! m_e      kg
      real(wp),parameter:: mproton    =  1.67262158e-27_wp    ! m_p      kg
      real(wp),parameter:: h_plank    =  6.62606876e-34_wp    ! h        j*s
      real(wp),parameter:: hbar_plank =  h_plank/(m_2pi)      ! h/2pi    j*s
      real(wp),parameter:: r_gas      =  8.314472_wp          ! r        j*mol(-1)*k(-1)
      real(wp),parameter:: k_boltz    =  1.3806503e-23_wp     ! k        j*k(-1)
      real(wp),parameter:: eV         =  eplus                ! ev       j
      real(wp),parameter:: faraday    =  96485.3415_wp        ! f        c*mol(-1)
      real(wp),parameter:: Angstrom1  =  1.0e-10_wp           !          m
      real(wp),parameter:: cal        =  4.1868_wp            !          j
!
!     Some derived constants used in MD code
      real(wp),parameter:: conv4   =  1000.0_wp*cal/(navo*k_boltz)
      real(wp),parameter:: conv4i  =  1.0_wp/conv4
      real(wp),parameter:: epstar  =  1000.0_wp*cal/(navo*k_boltz)
      real(wp),parameter,private:: ulength  =  Angstrom1
      real(wp),parameter:: p_epsilon  = k_boltz
      real(wp),parameter:: conv5  =  (ulength**2)/p_epsilon
END MODULE phys_const_mod

!!>INCLUDE 'global_vars_mod_KMC.f90'

MODULE global_vars_mod
   USE precision_mod
   USE phys_const_mod
   real(wp),parameter:: Angstrom = 1.0_wp  ! Angstroms/nm
   real(wp),parameter:: ulength = 1.0e-10_wp*Angstrom
   real(wp),parameter:: uenergy = eV
   real(wp),parameter:: upressure = uenergy/(ulength**3)
   real(wp),parameter:: torr2Pa = 101325.0_wp/760.0_wp
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = k_boltz/eV   ! 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp),parameter:: eVTokJmol = eV*navo/1000.0_wp
   real(wp),parameter:: kJmolToeV = 1.0_wp/eVTokJmol
   real(wp):: kBoltzT
   integer:: natom, natom_max, nfixed
   integer:: kmc_step, nrelax
END MODULE global_vars_mod

!!>INCLUDE 'files_mod.f90'

MODULE files_mod
   public get_free_file_unit, myflush,new_frame_file
CONTAINS

   SUBROUTINE get_free_file_unit(iu)
      implicit none
      character(len=*),parameter:: sub_name='get_free_file_unit'
      integer,intent(out):: iu
      integer,parameter:: istart=11
      integer,parameter:: unit_max=100000
      integer:: ios
      logical:: lopen,lexists
      do iu = istart,unit_max
         inquire(unit=iu, iostat=ios, exist=lexists, opened=lopen)
         if (ios == 0) then
            if (lexists .and. .not.lopen) RETURN
         else
            write(*,*) sub_name,': error, iostat = ',ios
            STOP
         end if
      end do
      write(*,*) sub_name,': error, no free units'
      STOP
   END SUBROUTINE get_free_file_unit


   SUBROUTINE new_frame_file(iu,base_name,iframe,file_type)
      integer,intent(in):: iframe
      integer,intent(out):: iu
      character(len=*),intent(in):: base_name
      character(len=3),intent(in):: file_type
      character(len=32):: ctmp
      character(len=132):: frame_file
      integer:: ios
      logical:: lexists,lopen
      call get_free_file_unit(iu)
      write(unit=ctmp,fmt='(i6.6)') iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.'//trim(file_type)
      inquire(file=trim(frame_file),iostat=ios,exist=lexists,opened=lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit=iu,file=trim(adjustl(frame_file)),form='formatted')
   END SUBROUTINE new_frame_file


   SUBROUTINE myflush(iu)
      integer,intent(in):: iu
      character(132):: fn
      inquire(unit=iu,name=fn)
      close(iu)
      open(unit=iu,file=trim(fn),position='append')
   END SUBROUTINE myflush

END MODULE files_mod

!!>INCLUDE 'sort_mod.f90'

MODULE sort_mod
   implicit none
   ! Some sorting routines modified from Numerical Recipes
   interface qsort
      module procedure qsort0,qsort1
   end interface qsort
CONTAINS


   SUBROUTINE qsort0(n,V)
      integer,intent(in):: n
      integer,intent(inout):: V(:)
! Sorts an array V(1:n) into ascending numerical order using the Quicksort algorithm.
! Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
! auxiliary storage.
      integer,parameter:: M = 7,NSTACK = 50
      integer:: i,ir,j,jstack,k,L,istack(NSTACK)
      integer:: a,temp
      jstack = 0
      L = 1
      ir = n
1     if ((ir - L) < M) then
         do j = L + 1,ir
            a = V(j)
            do i = j - 1,1, -1
               if (V(i) <= a) GOTO 2
               V(i + 1) = V(i)
            end do
            i = 0
2           V(i + 1) = a
         end do
         if (jstack == 0) RETURN
         ir = istack(jstack)
         L = istack(jstack - 1)
         jstack = jstack - 2
      else
         k = (L + ir)/2
         temp = V(k)
         V(k) = V(L + 1)
         V(L + 1) = temp
         if (V(L + 1) > V(ir)) then
            temp = V(L + 1)
            V(L + 1) = V(ir)
            V(ir) = temp
         end if
         if (V(L) > V(ir)) then
            temp = V(L)
            V(L) = V(ir)
            V(ir) = temp
         end if
         if (V(L + 1) > V(L)) then
            temp = V(L + 1)
            V(L + 1) = V(L)
            V(L) = temp
         end if
         i = L + 1
         j = ir
         a = V(L)
3        continue
         i = i + 1
         if (V(i) < a) GOTO 3
4        continue
         j = j - 1
         if (V(j) > a) GOTO 4
         if (j < i) GOTO 5
         temp = V(i)
         V(i) = V(j)
         V(j) = temp
         GOTO 3
5        V(L) = V(j)
         V(j) = a
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
      !  (C) Copr. 1986-92 Numerical Recipes Software #1-0zV'n26) B3.
   END SUBROUTINE qsort0


   SUBROUTINE qsort1(n,V0,V)
      integer,intent(in):: n,V0(:)
      integer,intent(out):: V(:)
      V(1:n) = V0(1:n)
      call qsort0(n,V)
   END SUBROUTINE qsort1


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

END MODULE sort_mod

!!>INCLUDE 'atom_types_mod.f90'

MODULE atom_types_mod
   USE precision_mod, only: wp
   implicit none
   integer,parameter:: idangling = 0
   integer,parameter:: iSilicon = 1
   integer,parameter:: iOxygen = 2
   integer,parameter:: iOxygenH = 3
   integer,parameter:: iHydrogen = 4
   integer,parameter:: iOw = 5
   integer,parameter:: iHw = 6
   integer,parameter:: iO_CO2 = 7
   integer,parameter:: iC_CO2 = 8
   integer,parameter:: iN_N2 = 9
   integer,parameter:: iO_O2 = 10
   integer,parameter:: ntyp = 10
   integer:: n_atomtype(ntyp)
   character(4),parameter:: atom_name(0:ntyp) = (/ ' __ ', &
      'Si  ','O   ','OH  ','H   ','Ow  ', &
      'Hw  ','O_CO','C_CO','N_N2','O_O2' /)
   integer,parameter:: ncmax(0:ntyp) = (/ 0, &
                       4, 2, 2, 1, 2, &
                       1, 2, 4, 3, 2 /)
   real(wp),parameter:: mSi = 28.06_wp
   real(wp),parameter:: mOx = 16.00_wp
   real(wp),parameter:: mOH = 17.00_wp
   real(wp),parameter:: mH =  1.00_wp
   real(wp),parameter:: mN = 14.00_wp
   real(wp),parameter:: mC = 12.00_wp
   real(wp),parameter:: amass(0:ntyp) = (/ 0.0_wp, &
                        mSi, mOx, mOH, mH, mOx, &
                         mH, mOx,  mC, mN, mOx /)
   integer,allocatable:: atom(:),atom2(:)
   logical,allocatable:: latom(:)
   integer,allocatable:: Ox_atom(:,:),N2_atom(:,:),CO2_atom(:,:),GAS_atom(:,:),h2o_atom(:,:)
   integer,private:: ido
   character(len=1), parameter,private:: upper_case(0:255) = (/ (achar(ido),   &
   ido = 0,96), (achar(ido - 32), ido = 97,122), (achar(ido), ido = 123,255) /)

   interface count_atoms
      module procedure count_atom_types, count_atom_types_list
   end interface count_atoms

CONTAINS

   PURE FUNCTION name2atom(str)
      integer:: name2atom
      character(*),intent(in):: str
      integer:: i,j
      outer: do i = 0,ntyp
         do j = 1,len_trim(str)
            if (upper_case(iachar( str(j:j) )) /= &
                upper_case(iachar(atom_name(i)(j:j))) ) then
               cycle outer
            end if
         end do
         name2atom = i
         RETURN
      end do outer
      if (i > ntyp) name2atom = -1
   END FUNCTION name2atom

   ELEMENTAL FUNCTION atomfn(i)
      integer:: atomfn
      integer,intent(in):: i
      if(i > 0) then
         atomfn = atom(i)
      else
         atomfn = 0
      end if
   END  FUNCTION atomfn

   SUBROUTINE count_atom_types(ifirst,ilast,atom,n_atom)
      integer,intent(in):: ifirst,ilast,atom(:)
      integer,intent(out):: n_atom(:)
      integer:: i,ia
      n_atom = 0
      do i = ifirst,ilast
         ia = atom(i)
         n_atom(ia) = n_atom(ia) + 1
      end do
   END SUBROUTINE count_atom_types

   SUBROUTINE count_atom_types_list(nl,lst,atom,n_atom)
      integer,intent(in):: nl,lst(:),atom(:)
      integer,intent(out):: n_atom(:)
      integer:: i,ia
      n_atom = 0
      do i = 1,nl
         ia = atom(lst(i))
         n_atom(ia) = n_atom(ia) + 1
      end do
   END SUBROUTINE count_atom_types_list

   SUBROUTINE write_atoms(iu,n_atom)
      integer,intent(in):: iu,n_atom(:)
      integer:: i
      do i = lbound(n_atom,1),ubound(n_atom,1)
         write(iu,'(i4,1x,a4,1x,i7)') i,atom_name(i),n_atom(i)
      end do
      write(iu,*)
   END SUBROUTINE write_atoms

END MODULE atom_types_mod

!!>INCLUDE 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod
   USE global_vars_mod, only: Angstrom
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
   integer,parameter:: NPER = 3
   logical,parameter:: box_periodic(ndim) = (/ ( (i <= nper), i=1,ndim ) /)
   real(wp):: boxl(NDIM),boxli(NDIM),boxl2(NDIM)
!
   TYPE CUBOID
      real(wp):: v1(3),v2(3)
   END TYPE CUBOID
!
   interface pbc
      module procedure  pbc_v, pbc_a, pbc_a2 , pbc_n_v, pbc_n_a, pbc_n_a2
   end interface
CONTAINS

   SUBROUTINE SET_BOX_SIZE(L)
      real(wp),intent(in):: L(NDIM)
      boxl = L
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
   END SUBROUTINE SET_BOX_SIZE

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

   PURE FUNCTION nearest_image(r,ro)
      real(wp):: nearest_image(NDIM)
      real(wp),intent(in):: r(:),ro(:)
      real(wp):: dr(3)
      dr = ro - r
      nearest_image = r
      where (abs(dr) > boxl2) nearest_image = r - sign(boxl,r)
   END FUNCTION nearest_image


   SUBROUTINE READ_PDB_HEADER(iu,natom,boxl)
      integer,intent(in):: iu
      integer,intent(inout):: natom
      real(wp),intent(inout)::boxl(:)
      character(6):: c6
      read(iu,*) c6,natom
      read(iu,*) c6,boxl
      boxl = boxl/Angstrom
   END SUBROUTINE READ_PDB_HEADER

   SUBROUTINE WRITE_PDB_HEADER(iu,natom,boxl)
      integer,intent(in):: iu,natom
      real(wp),intent(in):: boxl(:)
      write(iu,*) 'REMARK',natom
      write(iu,'(a6,1x,3(1x,f8.4))') 'REMARK',boxl*Angstrom
   END SUBROUTINE WRITE_PDB_HEADER


   SUBROUTINE WRITE_XYZ_HEADER(iu,natom,boxl,comment)
      integer,intent(in):: iu,natom
      real(wp),intent(in):: boxl(:)
      character(*),intent(in),optional:: comment
      write(iu,*) natom
      if (present(comment)) then
         write(iu,'(3(f0.8,1x),a)') boxl*Angstrom,trim(comment)
      else
         write(iu,'(3(f0.8,1x))') boxl*Angstrom
      end if
   END SUBROUTINE WRITE_XYZ_HEADER

   SUBROUTINE READ_XYZ_HEADER(iu,natom,boxl)
      integer,intent(in):: iu
      integer,intent(inout):: natom
      real(wp),intent(inout)::boxl(:)
      read(iu,*) natom
      read(iu,*) boxl
      boxl = boxl/Angstrom
   END SUBROUTINE READ_XYZ_HEADER

   SUBROUTINE WRITE_XYZ(iu,ifirst,ilast,atom,r_xyz)
      USE atom_types_mod, only: atom_name
      integer,intent(in):: iu,ifirst,ilast,atom(:)
      real(wp),intent(in):: r_xyz(:,:)
      integer:: i,j
      do i = ifirst,ilast
         write(iu,'(a4,3(1x,f14.8))') atom_name(atom(i)),(r_xyz(i,j)*Angstrom,j = 1,NDIM)
      end do
      !write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu,natom,atom,r_xyz)
      USE atom_types_mod, only: name2atom
      integer,intent(in):: iu,natom
      integer,intent(inout):: atom(:)
      real(wp),intent(inout):: r_xyz(:,:)
      integer:: i
      real(wp):: r(NDIM)
      character(32):: ctmp
      do i = 1,natom
         read(iu,*) ctmp,r
         atom(i) = name2atom(trim(ctmp))
         r_xyz(i,1:NDIM) = r/Angstrom
      end do
   END SUBROUTINE READ_XYZ

   SUBROUTINE WRITE_XYZ_FILE(iu,ifirst,ilast,atom,r_xyz,boxl,comment,box_corners)
      integer,intent(in):: iu,ifirst,ilast,atom(:)
      real(wp),intent(in):: r_xyz(:,:),boxl(:)
      character(*),intent(in),optional:: comment
      type(cuboid),intent(in),optional:: box_corners
      integer:: nat
      logical:: print_box
      print_box = .false.
      if (present(box_corners)) print_box = .true.
      nat = ilast-ifirst+1
      if (print_box) nat = nat + 8
      if(present(comment))then
         call WRITE_XYZ_HEADER(iu,nat,boxl,comment)
      else
         call WRITE_XYZ_HEADER(iu,nat,boxl)
      end if
      call WRITE_XYZ(iu,ifirst,ilast,atom,r_xyz)
      if (print_box) call WRITE_BOX_XYZ(iu,box_corners%v1*Angstrom,box_corners%v2*Angstrom,'Cr')
      write(iu,*)
      close(iu)
   END SUBROUTINE WRITE_XYZ_FILE

   SUBROUTINE write_xyz_sphere(iu,r_xyz,atoms,ifirst,ilast,ro,rad)
      USE atom_types_mod,only: atom_name
      integer,intent(in):: iu,ifirst,ilast,atoms(:)
      real(wp),intent(in):: r_xyz(:,:),ro(:),rad
      real(wp):: r2,dr(3),rv(3)
      integer:: i
      r2 = rad*rad
      do i = ifirst,ilast
         dr = ro - r_xyz(i,1:3)
         call pbc(dr)
         if (dot_product(dr,dr) < r2) then
            dr = ro - r_xyz(i,1:3)
            rv = r_xyz(i,:)
            where (abs(dr) > boxl2) rv = r_xyz(i,:) - sign(boxl,r_xyz(i,:))
            write(iu,'(a,3(f12.6))') atom_name(atoms(i)),rv*Angstrom
         end if
      end do
   END SUBROUTINE write_xyz_sphere


   SUBROUTINE WRITE_PDB(iu,ifirst,ilast,atom,r_xyz)
      USE atom_types_mod, only: atom_name
      integer,intent(in):: iu,ifirst,ilast,atom(:)
      real(wp),intent(in):: r_xyz(:,:)
      integer:: i
      do i = ifirst,ilast
         write(iu,110)'HETATM',i,atom_name(atom(i)),'   ',i,r_xyz(i,1:NDIM)*Angstrom
      end do
110   format(a6,i5,1x,a4,1x,a3,i6,4x,3f8.3)
   END SUBROUTINE WRITE_PDB

   SUBROUTINE READ_PDB(iu,natom,atom,r_xyz)
      USE atom_types_mod, only: name2atom
      integer,intent(in):: iu,natom
      integer,intent(inout):: atom(:)
      real(wp),intent(inout):: r_xyz(:,:)
      integer:: i,itmp
      real(wp):: rv(NDIM)
      character(32):: ctmp,c6*6
      do i = 1,natom
         read (iu,'(a6)',advance = 'no') c6
         read (iu,*) itmp,ctmp,itmp,rv
         atom(i) = name2atom(trim(ctmp))
         r_xyz(i,1:NDIM) = rv/Angstrom
      end do
   END SUBROUTINE READ_PDB

   SUBROUTINE READ_PDB_MOL(iu,N,r_xyz,atom)
      integer,intent(in):: iu,N
      real(wp),intent(out):: r_xyz(:,:)
      integer,intent(out):: atom(:)
      integer:: itmp
      character(6):: ctmp1
      read(iu,*) ctmp1, itmp
      if (itmp /= N) STOP 'incorrect no. of atoms in file'
      read(iu,*) ctmp1
      call READ_PDB(iu,N,atom,r_xyz)
      close(iu)
   END SUBROUTINE READ_PDB_MOL

   SUBROUTINE WRITE_BOX_PDB(iu,rBot,rTop,ia,ch)
      integer,intent(in):: iu,ia
      real(wp),intent(in):: rBot(3),rTop(3)
      character(*):: ch
      write(iu,110)'HETATM',ia+0,trim(ch),'   ',ia+0, rBot(1), rBot(2), rBot(3)
      write(iu,110)'HETATM',ia+1,trim(ch),'   ',ia+1, rTop(1), rBot(2), rBot(3)
      write(iu,110)'HETATM',ia+2,trim(ch),'   ',ia+2, rBot(1), rTop(2), rBot(3)
      write(iu,110)'HETATM',ia+3,trim(ch),'   ',ia+3, rTop(1), rTop(2), rBot(3)
      write(iu,110)'HETATM',ia+4,trim(ch),'   ',ia+4, rBot(1), rBot(2), rTop(3)
      write(iu,110)'HETATM',ia+5,trim(ch),'   ',ia+5, rTop(1), rBot(2), rTop(3)
      write(iu,110)'HETATM',ia+6,trim(ch),'   ',ia+6, rBot(1), rTop(2), rTop(3)
      write(iu,110)'HETATM',ia+7,trim(ch),'   ',ia+7, rTop(1), rTop(2), rTop(3)
      write(iu,"(a6,2i5)")'CONECT',ia+0,ia+1
      write(iu,"(a6,2i5)")'CONECT',ia+0,ia+2
      write(iu,"(a6,2i5)")'CONECT',ia+0,ia+4
      write(iu,"(a6,2i5)")'CONECT',ia+1,ia+3
      write(iu,"(a6,2i5)")'CONECT',ia+1,ia+5
      write(iu,"(a6,2i5)")'CONECT',ia+2,ia+3
      write(iu,"(a6,2i5)")'CONECT',ia+2,ia+6
      write(iu,"(a6,2i5)")'CONECT',ia+3,ia+7
      write(iu,"(a6,2i5)")'CONECT',ia+4,ia+5
      write(iu,"(a6,2i5)")'CONECT',ia+4,ia+6
      write(iu,"(a6,2i5)")'CONECT',ia+5,ia+7
      write(iu,"(a6,2i5)")'CONECT',ia+6,ia+7
110   format(a6,i5,1x,a4,1x,a3,i6,4x,3f8.3)
   END SUBROUTINE WRITE_BOX_PDB

   SUBROUTINE WRITE_BOX_XYZ(iu,rBot,rTop,ch)
      integer,intent(in):: iu
      real(wp),intent(in):: rBot(3),rTop(3)
      character(*):: ch
      write(iu,'(a,3f12.6)') trim(ch), rBot(1), rBot(2), rBot(3)
      write(iu,'(a,3f12.6)') trim(ch), rTop(1), rBot(2), rBot(3)
      write(iu,'(a,3f12.6)') trim(ch), rBot(1), rTop(2), rBot(3)
      write(iu,'(a,3f12.6)') trim(ch), rTop(1), rTop(2), rBot(3)
      write(iu,'(a,3f12.6)') trim(ch), rBot(1), rBot(2), rTop(3)
      write(iu,'(a,3f12.6)') trim(ch), rTop(1), rBot(2), rTop(3)
      write(iu,'(a,3f12.6)') trim(ch), rBot(1), rTop(2), rTop(3)
      write(iu,'(a,3f12.6)') trim(ch), rTop(1), rTop(2), rTop(3)
   END SUBROUTINE WRITE_BOX_XYZ

END MODULE coordinates_mod

!!>INCLUDE 'connectivity_mod.f90'

MODULE connectivity_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom, nfixed
   USE atom_types_mod
   implicit none
   integer,allocatable:: proximity(:,:)
!
   interface change_connectivity
      module procedure  set_proximity, set_proximity2
   end interface change_connectivity
!
CONTAINS

   SUBROUTINE write_atom_conn(iu,ifirst,ilast,atom,conect)
      USE atom_types_mod,only: ncmax,name2atom
!     write topology/connectivity information to file unit iu
      integer,intent(in):: iu,ifirst,ilast,atom(:),conect(:,:)
      integer:: i,j,k
      write(iu,*) ilast-ifirst+1
      do i = ifirst,ilast
         k = count(conect(i,:) > 0)
         write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
            i,atom_name(atom(i)),k,(conect(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE write_atom_conn


   SUBROUTINE read_atom_conn(iu,natom,atom,conect)
      USE atom_types_mod,only: ncmax,name2atom
!     Read topology/connectivity information from file unit iu
      integer,intent(in):: iu
      integer,intent(inout):: natom,atom(:),conect(:,:)
      integer:: i,j,ic,ii
      character(32):: ctmp
      read(iu,*) natom
      do i = 1,natom
         read(iu,*) ii,ctmp,ic,(conect(i,j),j = 1,ncmax(name2atom(ctmp)))
         atom(i) = name2atom(ctmp)
         if ( atom(i) < 1 ) STOP 'error reading atom type'
      end do
      rewind(iu)
   END SUBROUTINE read_atom_conn


   SUBROUTINE write_conect_pdb(iu,ifirst,ilast,conect)
      ! write CONECT records to a PDB file
      use atom_types_mod
      integer,intent(in):: iu,ifirst,ilast,conect(:,:)
      integer:: i,j
      do i = ifirst,ilast
         if ( all(conect(i,:) == 0) ) cycle
         write(iu,"(a6,9i5)")'CONECT',i,(conect(i,j),j = 1,ncmax(atom(i)))
      end do
      call flush(iu)
   END SUBROUTINE write_conect_pdb


   SUBROUTINE read_conect_pdb(iu,natom,conect)
      ! read natom CONECT records from a PDB file
      use atom_types_mod
      integer,intent(in):: iu,natom
      integer,intent(inout):: conect(:,:)
      integer:: i,ia,j
      character(len=80):: line,c5*5
      do ia = 1,natom
         read(iu,'(a80)') line
         c5 = line(7:11)
         read(unit=c5,fmt=*) i
         do j = 1,ncmax(atom(i))
            c5 = line(6+5*j+1:6+5*j+5)
            read(unit=c5,fmt=*) conect(i,j)
         end do
      end do
   END SUBROUTINE read_conect_pdb


   SUBROUTINE write_conn(iu,ifirst,ilast,conect,line)
      USE atom_types_mod
      integer,intent(in):: iu,ifirst,ilast,conect(:,:)
      character(len=*),intent(in):: line
      integer:: j,ia
      ! write out connectivity
      write(iu,*) ilast-ifirst+1
      write(iu,*) line
      do ia = ifirst,ilast
         write(iu,'(9(1x,i7))') ia,(conect(ia,j),j=1,ncmax(atom(ia)))
      end do
      call flush(iu)
   END SUBROUTINE write_conn


   SUBROUTINE read_conn(iu,conect)
      USE atom_types_mod
      integer,intent(in):: iu
      integer,intent(inout):: conect(:,:)
      integer:: nat,i,j,ia
      character(len=132):: line
      ! read in connectivity
      read(iu,*) nat
      read(iu,*) line
      do i = 1,nat
         read(iu,*) ia,(conect(ia,j),j=1,ncmax(atom(ia)))
      end do
   END SUBROUTINE read_conn


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


   SUBROUTINE set_proximity2(iat,iold,inew,proximity)
!     Where the value of proximity is iold set it to inew.
!     Used to break/switch bonds.
!     e.g. the following move requires two call to set_proximity
!
!     A--B         A  B      set_proximity(A,B,C)  switch A-B to A-C
!            -->    \        set_proximity(B,A,0)  break B-A
!        C           C
!
      integer,intent(in):: iat,iold,inew
      integer,intent(inout):: proximity(:,:)
      integer:: i
      do i = 1,ncmax(atom(iat))
         if (proximity(iat,i) == iold) then
            proximity(iat,i) = inew
            exit
         end if
         if (i == ncmax(atom(iat))) then
            print *,'>> ',iat,iold,inew
            STOP 'set_proximity: ERROR'
         end if
      end do
   END SUBROUTINE set_proximity2


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
   END FUNCTION is_dangling_group


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
         if (proximity(iSi2,i) == iO1) cycle
         iat = proximity(iSi2,i)
         if (ANY(proximity(iat,:) == iSi1)) then
            is_2Si_ring = .TRUE.
            RETURN
         end if
      end do
   END FUNCTION is_2Si_ring


    FUNCTION OH_groups_OK(iOH1,iOH2)
      logical:: OH_groups_OK
      integer,intent(in):: iOH1,iOH2
      integer:: iSi1,iSi2
      OH_groups_OK = .FALSE.
      iSi1 = proximity(iOH1,1)
!print *,'OH_groups_OK: ',iOH1,iOH2
!print *,'iSi1 = ',iSi1
      if (iSi1 == 0) iSi1 = proximity(iOH1,2)
!print *,'atom(iSi1) = ',atom(iSi1)
      !if (atom(iSi1) /= iSilicon) iSi1 = proximity(iOH1,2)
      iSi2 = proximity(iOH2,1)
      if (iSi2 == 0) iSi2 = proximity(iOH2,2)
      !if (atom(iSi2) /= iSilicon) iSi2 = proximity(iOH2,2)
      if (iSi1 == iSi2) RETURN
if (iSi1 <= nfixed .or. iSi2 <= nfixed) then
   OH_groups_OK = .TRUE.
   RETURN
end if
      if (.NOT.nearest_neighbor2(iSi1,iSi2)) OH_groups_OK = .TRUE.
   END FUNCTION OH_groups_OK


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
   END FUNCTION noh


   SUBROUTINE check_proximity(ifirst,ilast,consistent,ierr)
!     Check that if atom i is bonded to k then
!     atom k is bonded to i, if everything is ok
!     then consistent is .true.
      integer,intent(in):: ifirst,ilast
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (k == 0) cycle
            if ( count(proximity(k,:) == i) /= 1) then
               consistent = .false.
               ierr = i
               RETURN
            end if
         end do
      end do
   END SUBROUTINE check_proximity


   SUBROUTINE delete_atom(i)
!     Delete atom i from the system and renumber the
!     last atom (natom) to number i
      USE coordinates_mod
      integer,intent(in):: i
      integer:: j,k
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         where( proximity(k,:) == i ) proximity(k,:) = 0
      end do
      rxyz(i,1:NDIM) = rxyz(natom,1:NDIM)
      proximity(i,:) = proximity(natom,:)
      atom(i) = atom(natom)
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         where( proximity(k,:) == natom ) proximity(k,:) = i
      end do
      rxyz(natom,1:NDIM) = 0.0_wp
      proximity(natom,:) = 0
      atom(natom) = 0
      natom = natom - 1
   END SUBROUTINE delete_atom


   SUBROUTINE delete_atom_list(n,ind)
      USE sort_mod, only: qsort
      integer,intent(in):: n
      integer,intent(inout):: ind(:)
      integer:: i
      call qsort(n,ind)
      do i = n,1,-1
         call delete_atom(ind(i))
      end do
   END SUBROUTINE delete_atom_list


   SUBROUTINE delete_atom_list2(n,ind,keep)
      USE sort_mod, only: qsort
      integer,intent(in):: n
      integer,intent(inout):: ind(:),keep(:)
      integer:: i,j
      call qsort(n,ind)
      do i = n,1,-1
         ! where (keep == natom) keep = ind(i)
         do j = 1,size(keep)
            if (keep(j) == natom) then
               keep(j) = ind(i)
               exit
            end if
         end do
         call delete_atom(ind(i))
      end do
   END SUBROUTINE delete_atom_list2


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
            STOP 'not a dangling group'
         end if
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE delete_group0

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
            STOP 'not a dangling group'
         end if
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE delete_group

   SUBROUTINE delete_group2(iSi)
!     Delete atom iSi and its nearest neighbors up to order 2
!     e.g. delete Si(OH)4 cluster
      USE sort_mod,only:shell
      integer,intent(in):: iSi
      integer:: k,j,i,m
      integer:: ilst(9),ni
      ni = 0
      ilst = 0
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == 0) cycle
         if (atom(j) /= iOxygenH) STOP 'not a dangling group'
         do m = 1,ncmax(atom(j))
            i = proximity(j,m)
            if (i == iSi) cycle
            if (i == 0) cycle
            !if (atom(i) /= iHydrogen) STOP 'ERROR: delete_group2' !only for debugging
            ni = ni + 1
            ilst(ni) = i
         end do
         ni = ni + 1
         ilst(ni) = j
      end do
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE delete_group2


   SUBROUTINE delete_group_C2H5(iSi,iOc)
!     Delete C2H5 cluster
      USE sort_mod,only:shell
      USE atom_types_mod
      integer,intent(inout):: iSi,iOc
      integer:: k,j,i,ii,jj,kk
      integer:: ilst(7),ni
      ni = 0
      ilst = 0
      loop_1: do i = 1,ncmax(atom(iOc))
         ii = proximity(iOc,i)
         if (ii == iSi) cycle loop_1
         if (ii == 0) STOP 'error: delete_group_C2H5'
         ni = ni + 1
         ilst(ni) = ii
         loop_2: do j = 1,ncmax(atom(ii))
            jj = proximity(ii,j)
            if (jj == iOc) cycle loop_2
            ni = ni + 1
            ilst(ni) = jj
            loop_3: do k = 1,ncmax(atom(jj))
               kk = proximity(jj,k)
               if (kk == ii) cycle loop_3
               ni = ni + 1
               ilst(ni) = kk
            end do loop_3
         end do loop_2
      end do loop_1
      call shell(ni,ilst)
      do i = ni,1,-1
         if (iOc == natom) then ! update the oxygen number
            iOc = ilst(i)
         else if (iSi == natom) then ! update the silicon number
            iSi = ilst(i)
         end if
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE delete_group_C2H5


   SUBROUTINE delete_group_OH(iSi,ioh)
      ! delete only OH group attached with Si
      USE sort_mod, only:shell
      USE global_vars_mod
      USE atom_types_mod
      integer,intent(inout):: iSi
      integer,intent(in):: ioh
      integer:: i,ii,ilst(2),ni
      ni = 1
      ilst(ni) = ioh
      do i = 1,ncmax(atom(ioh))
         ii = proximity(ioh,i)
         if (ii == iSi) cycle
         ni = ni + 1
         ilst(ni) = ii
      end do
      call shell(ni,ilst)
      do i = ni,1,-1
         if (iSi == natom) then ! update the silicon number
            iSi = ilst(i)
         end if
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE delete_group_OH


   SUBROUTINE delete_Hydrogen(iSi,iOH)
      USE sort_mod, only:shell
      implicit none
      integer,intent(inout):: iSi,iOH
      integer:: i,ii
      ! delete only H  attached with OH
      do i = 1,ncmax(atom(iOH))
         ii = proximity(iOH,i)
         if (ii == 0) cycle
         if (atom(ii) == iHydrogen) exit
      end do
      if (i > ncmax(atom(iOH))) STOP 'delete_Hydrogen: error'
      if (iOH == natom) then ! update the oxygen number
         iOH = ii
      else if (iSi == natom) then ! update the silicon number
         iSi = ii
      end if
      call delete_atom(ii)
   END SUBROUTINE delete_Hydrogen


   SUBROUTINE find_1st_atomtype(ia,itype,ib)
      integer,intent(in):: ia, itype
      integer,intent(out):: ib
      integer:: i,j
! find ib, the 1st atom of type itype connected to ia
!
!   |
! --ia--ib (= itype)
!   |
!
      ib = 0
      do i = 1, ncmax(atom(ia))
         j = proximity (ia,i)
         if ( j == 0 ) cycle
         if (atom(j) == itype) then
            ib = j
            return
         end if
      end do
      stop 'error: atom not found'
   END SUBROUTINE find_1st_atomtype


   SUBROUTINE list_connected_atoms1(iSi,itype,n,Lst)
      integer,intent(in):: iSi, itype
      integer,intent(out):: n, Lst(:)
      integer:: i,j
!
!    C2H5 C2H5
!      |  |
!      O  O  O-C2H5
!       \ | /
!         Si
!         |
!
      n = 0
      Lst = 0
      do i = 1, ncmax(atom(iSi))
         j = proximity (iSi,i)
         if ( j == 0 ) cycle
         if (atom(j) == itype) then
            n = n + 1
            Lst(n) = j
         end if
      end do
    END SUBROUTINE list_connected_atoms1


    SUBROUTINE list_connected_atoms2(iSi,itype1,itype2,n1,Lst1,n2,Lst2)
!
!      H C2H5                H  H
!      |  |                  |  |
!      O  O  O-C2H5          O  O  O-C2H5
!       \ | /                 \ | /
!         Si                    Si
!         |                     |
!
      integer,intent(in):: iSi, itype1,itype2
      integer,intent(out)::n1,Lst1(:),n2,Lst2(:)
      integer:: i,j
      n1 = 0
      n2 = 0
      Lst1 = 0
      Lst2 = 0
      do i = 1, ncmax(atom(iSi))
         j = proximity(iSi, i)
         if (atom(j) == itype1) then
            n1 = n1 + 1
            Lst1(n1) = j
         else if (atom(j) == itype2) then
            n2 = n2 + 1
            Lst2(n2) = j
         end if
      end do
   END SUBROUTINE list_connected_atoms2


   SUBROUTINE add_nn_list(ia,list,nl)
      integer,intent(in):: ia
      integer,intent(inout):: list(:),nl
      integer:: j,jj
      do jj = 1, ncmax(atom(ia))
         j = proximity(ia,jj)
         if (j == 0) cycle
         nl = nl + 1
         list(nl) = j
      end do
   END SUBROUTINE add_nn_list


   SUBROUTINE add_nn_list_uniq(ia,list,nl)
      integer,intent(in):: ia
      integer,intent(inout):: list(:),nl
      integer:: j,jj
      do jj = 1, ncmax(atom(ia))
         j = proximity(ia,jj)
         if (j == 0) cycle
         if (in_list(j,list,nl)) cycle
         nl = nl + 1
         list(nl) = j
      end do
   contains
      LOGICAL PURE FUNCTION in_list(k,list,nl)
         integer,intent(in):: k,list(:),nl
         integer:: i
         in_list = .FALSE.
         do i = 1,nl
            if (list(i) == k) then
               in_list = .TRUE.
               RETURN
            end if
         end do
      END FUNCTION in_list
   END SUBROUTINE add_nn_list_uniq


   SUBROUTINE swap_atom(a,b)
      USE coordinates_mod
!      USE connectivity_mod
      integer,intent(in):: a,b
      integer:: itmp,proxtmp(size(proximity,2))
      integer:: list(2*size(proximity,2) + 2),nl,i,jj,j,il
      real(wp):: rtmp(NDIM)
      nl = 2
      list(1) = a
      list(2) = b
      call add_nn_list_uniq(a,list,nl)
      call add_nn_list_uniq(b,list,nl)
      proxtmp = proximity(a,:)
      proximity(a,:) = proximity(b,:)
      proximity(b,:) = proxtmp
      do il = 1,nl
         i = list(il)
         do jj = 1,ncmax(atom(i))
            j = proximity(i,jj)
            if (j == 0) cycle
            if (j == a) then
               proximity(i,jj) = b
            else if (j == b) then
               proximity(i,jj) = a
            end if
         end do
      end do
      rtmp = rxyz(a,:)
      rxyz(a,:) = rxyz(b,:)
      rxyz(b,:) = rtmp
      itmp = atom(a)
      atom(a) = atom(b)
      atom(b) = itmp
   END SUBROUTINE swap_atom


   SUBROUTINE swap_atoms(a,b)
      USE coordinates_mod, only: rxyz
      integer,intent(in):: a,b
      integer:: j,k,atomtemp,proxtemp(4)
      real(wp):: rxyztemp(3)
      do j = 1,ncmax(atom(a))
         k = proximity(a,j)
         if (k == 0) cycle
         where( proximity(k,:) == a ) proximity(k,:) = b
      end do
      rxyztemp = rxyz(a,1:3)
      rxyz(a,1:3) = rxyz(b,1:3)
      rxyz(b,1:3) = rxyztemp

      proxtemp = proximity(a,:)
      proximity(a,:) = proximity(b,:)
      proximity(b,:) = proxtemp

      atomtemp = atom(a)
      atom(a) = atom(b)
      atom(b) = atomtemp

      do j = 1,ncmax(atom(a))
         k = proximity(a,j)
         if (k == 0) cycle
         where( proximity(k,:) == b ) proximity(k,:) = a
      end do
   END SUBROUTINE swap_atoms


   SUBROUTINE write_atom_info(iu,ia)
      integer,intent(in):: iu,ia
      integer:: j,jj
      write(iu,'(i6,2x,a2,2x,4(1x,i6),4x,4(1x,a4))')ia,atom_name(atom(ia)),proximity(ia,:),atom_name(atomfn(proximity(ia,:)))
      do jj = 1,ncmax(atom(ia))
         j = proximity(ia,jj)
         if(j==0)cycle
         write(iu,'(i6,2x,a2,2x,4(1x,i6),4x,4(1x,a4))')j,atom_name(atom(j)),proximity(j,:),atom_name(atomfn(proximity(j,:)))
      end do
   END SUBROUTINE write_atom_info

END MODULE connectivity_mod

!!>INCLUDE 'nlist_mod3.f90'

MODULE NLIST_MOD
   USE precision_mod
   USE coordinates_mod
   implicit none
   PRIVATE
   integer,PUBLIC:: ncelx,ncely,ncelz,ncelt
   real(wp),PUBLIC:: del_nl_x,del_nl_y,del_nl_z
   real(wp),private:: delxi,delyi,delzi,delmin
   real(wp),private:: bottom_corner(3)
   integer,allocatable,PUBLIC:: hoc(:),hoc_old(:)
   integer,allocatable,PUBLIC:: ll_old(:),ll(:),lr(:),lr_old(:)
   PUBLIC:: NEIGCELL,cell,INIT_NLIST,NEW_NLIST,print_cell,cell2xyz, icell, delete_nlist
   PUBLIC:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll,store_nlist,restore_nlist
!
!  ll(i)      : linked list particle i
!  hoc(ic)    : head of chain cell ic
!  ncelx,y,z  : number of cells in x, y or z direction
!  ncelt      : total number of cells
!
   interface cell2xyz
      module procedure icell2xyz, cell2ir
   end interface cell2xyz

   interface icell
      module procedure cell_ixyz, cell_ir
   end interface  icell

CONTAINS

   PURE FUNCTION nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   END FUNCTION

   SUBROUTINE INIT_NLIST(rlower,rupper,Rc,natom_max)
!
!          ________(b)
!         /|      /|
!        / |     / |  Neighbor list box
!       /__|____/  |
!       |  |____|__|
!       | /     |  /
!       |/      | /   (a) is position rlower = bottom_corner
!       |_______|/    (b) is position rupper
!      (a)
!
      real(wp),intent(in):: rlower(3),rupper(3),Rc
      integer,intent(in):: natom_max
      real(wp):: side(3)
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max),lr_old(natom_max))
      bottom_corner = rlower
      side = rupper - rlower
      ncelx = int(side(1)/Rc)
      ncely = int(side(2)/Rc)
      ncelz = int(side(3)/Rc)
      ncelt = ncelx*ncely*ncelz
      del_nl_x = side(1)/ncelx
      del_nl_y = side(2)/ncely
      del_nl_z = side(3)/ncelz
      delmin = min(del_nl_x,del_nl_y,del_nl_z)
      delxi = 1.0_wp/del_nl_x
      delyi = 1.0_wp/del_nl_y
      delzi = 1.0_wp/del_nl_z
      write(*,'(/a/a,3i6)') 'INIT_NLIST','ncel = ',ncelx,ncely,ncelz
      write(*,'(a,i9)') 'ncelt = ',ncelt
      write(*,'(a,3f12.6)') 'del = ',del_nl_x,del_nl_y,del_nl_z
      write(*,'(a,f12.6/)') 'delmin = ',delmin
      allocate(hoc(ncelt),hoc_old(ncelt))
      LL = 0
      HOC = 0
      LR = 0
   END SUBROUTINE INIT_NLIST


   SUBROUTINE DELETE_NLIST()
      if (allocated(ll_old)) deallocate(ll_old)
      if (allocated(ll)) deallocate(ll)
      if (allocated(lr)) deallocate(lr)
      if (allocated(lr_old)) deallocate(lr_old)
      if (allocated(hoc)) deallocate(hoc)
      if (allocated(hoc_old)) deallocate(hoc_old)
   END SUBROUTINE DELETE_NLIST


   PURE FUNCTION CELL(XYZ)
!     determines cell number for position x,y,z
      real(wp),intent(in)::  XYZ(:)
      integer:: CELL
      CELL = INT((XYZ(1) - bottom_corner(1))*delxi) &
           + INT((XYZ(2) - bottom_corner(2))*delyi)*ncelx &
           + INT((XYZ(3) - bottom_corner(3))*delzi)*ncelx*ncely + 1
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
         if (ic < 1 .or. ic > ncelt) STOP 'NEW_NLIST: cell number out of bounds'
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

   PURE SUBROUTINE icell2xyz(ic,ix,iy,iz)
      integer,intent(in):: ic
      integer,intent(out):: ix,iy,iz
      integer:: n2
      n2 = ncelx*ncely
      iz = ic/n2 + 1
      if (mod(ic,n2) == 0) iz = iz - 1
      iy = (ic - (iz - 1)*n2)/ncelx + 1
      if (mod(ic - (iz - 1)*n2,ncelx) == 0) iy = iy - 1
      ix = ic - (iy - 1)*ncelx - (iz - 1)*n2
   END SUBROUTINE icell2xyz

   PURE SUBROUTINE cell2ir(ic,ir)
      integer,intent(in):: ic
      integer,intent(out):: ir(3)
      integer:: n2
      n2 = ncelx*ncely
      ir(3) = ic/n2 + 1
      if (mod(ic,n2) == 0) ir(3) = ir(3) - 1
      ir(2) = (ic - (ir(3) - 1)*n2)/ncelx + 1
      if (mod(ic - (ir(3) - 1)*n2,ncelx) == 0) ir(2) = ir(2) - 1
      ir(1) = ic - (ir(2) - 1)*ncelx - (ir(3) - 1)*n2
   END SUBROUTINE cell2ir

   PURE FUNCTION cell_ixyz(ix,iy,iz)
      integer:: cell_ixyz
      integer,intent(in):: ix,iy,iz
      cell_ixyz = (iz - 1)*ncelx*ncely + (iy - 1)*ncelx + ix
   END FUNCTION cell_ixyz

   PURE FUNCTION cell_ir(ir)
      integer:: cell_ir
      integer,intent(in):: ir(3)
      cell_ir = (ir(3) - 1)*ncelx*ncely + (ir(2) - 1)*ncelx + ir(1)
   END FUNCTION cell_ir


   PURE FUNCTION BAD_CELL_NUM(ic)
      integer,intent(in):: ic
      logical:: BAD_CELL_NUM
      BAD_CELL_NUM = (ic < 1).or.(ic > ncelt)
   END FUNCTION BAD_CELL_NUM

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

!!>INCLUDE 'HKNonLattice2.f90'

MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster=0,n_cluster_old=0
   integer,allocatable:: atomL(:),atomL_old(:)

   interface Init_HKNonLattice
      module procedure Init_HKNonLattice_default
      module procedure Init_HKNonLattice_general
   end interface Init_HKNonLattice

CONTAINS

   SUBROUTINE Init_HKNonLattice_default(natom_max)
      integer,intent(in):: natom_max
      if (.not.allocated(atomL)) then
         allocate(atomL(natom_max),atomL_old(natom_max))
      else
         if (size(atomL) < natom_max) then
            deallocate(atomL,atomL_old)
            allocate(atomL(natom_max),atomL_old(natom_max))
         end if
      endif
      n_cluster = 0
      n_cluster_old = 0
      atomL = 0
      atomL_old = 0
   END SUBROUTINE Init_HKNonLattice_default


   SUBROUTINE Init_HKNonLattice_general(NumberOfNodes,NodeL)
      integer,intent(in):: NumberOfNodes
      integer,allocatable,intent(inout):: NodeL(:)
      if (.not.allocated(NodeL)) then
         allocate(NodeL(NumberOfNodes))
      end if
      NodeL = 0
   END SUBROUTINE Init_HKNonLattice_general


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
      nnmax = size(NodeNext,dim=2)
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
               NodeNextL(iNodeNextL) = NodeL(j)
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
      where (NodeLP1(2:iNodeLP) > NodeLP1(1:iNodeLP-1)) RelabLB(2:iNodeLP) = 1
      RelabL1(1:iNodeLP-1) = NodeLP1(2:iNodeLP)*RelabLB(2:iNodeLP)
      RelabL = 0
      RelabL(1) = NodeLP1(1)
      ii = 1
      do i = 1,iNodeLP-1
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
   END SUBROUTINE HKNonLattice


   SUBROUTINE cluster_count(n_node,n_cluster,NodeL,MaxCluster,ClusterCount)
      integer,intent(in):: n_node,n_cluster,NodeL(:)
      integer,intent(out):: maxcluster
      integer,allocatable,intent(out):: ClusterCount(:)
      integer:: i,L
      allocate(ClusterCount(n_cluster))
      ClusterCount = 0
      do i = 1,n_node
         L = NodeL(i)
         ClusterCount(L) = ClusterCount(L) + 1
      end do
      MaxCluster = maxloc(ClusterCount,1)
      ! some consistency checks
      do i = 1,n_cluster
         if (ClusterCount(i) /= count(NodeL(1:n_node) == i)) then
            stop 'error in cluster_count'
         end if
      end do
      if (sum(ClusterCount(1:n_cluster)) /= n_node) then
         stop 'error in cluster_count: wrong count'
      end if
   END SUBROUTINE cluster_count


! some routines used mainly for debugging
   SUBROUTINE analyse_cluster(n_node,NodeL,n_cluster,clusterC,cluster)
      integer,intent(in):: n_node,NodeL(:),n_cluster
      integer,intent(out):: clusterC(:),cluster(:,:)
      integer:: i,L
      clusterC = 0
      cluster = 0
      do i = 1,n_node
         L = NodeL(i)
         clusterC(L) = clusterC(L) + 1
         cluster(L,clusterC(L)) = i
      end do
      ! some consistency checks
      do i = 1,n_cluster
         if (clusterC(i) /= count(NodeL(1:n_node) == i)) then
            stop 'error in analyse_cluster'
         end if
      end do
      if (sum(clusterC(1:n_cluster)) /= n_node) then
         stop 'error in analyse_cluster: wrong count'
      end if
   END SUBROUTINE analyse_cluster

   SUBROUTINE print_cluster(iu,n_cluster,clusterC,cluster,n_node,NodeL)
      integer,intent(in):: iu,n_cluster,clusterC(:),cluster(:,:)
      integer,intent(in),optional:: n_node,NodeL(:)
      integer:: i
      do i = 1,n_cluster
         write(iu,*)'cluster = ',i,'count = ',clusterC(i)
         write(iu,'(10i7)') cluster(i,1:clusterC(i))
         if (present(n_node).and.present(NodeL)) then
            if (clusterC(i) /= count(NodeL(1:n_node) == i)) stop 'error in print_cluster'
         end if
      end do
   END SUBROUTINE print_cluster

END MODULE HKNonLattice_mod

!!>INCLUDE 'check_structure_mod_Silica.f90'

MODULE check_structure_mod
   USE global_vars_mod
   USE atom_types_mod
CONTAINS

   SUBROUTINE check_proximity2(ifirst,ilast,atoms,conect,consistent,ierr)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: atoms(:),conect(:,:)
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k,ia,a1,a2,i1,i2
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         ia = atoms(i)
         select case(ia)
         case(iSilicon)  ! Check Si is bonded to 4 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100   ! but not Si
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            if (a1 /= iSilicon) GOTO 100
            i2 = conect(i,2)  !  must be bonded to  2 Silicons
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (a2 /= iSilicon) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
         case(iOxygenH)  ! OH Oxygen must be bonded to Si
            i1 = conect(i,1)
            ! if (i1 == 0) GOTO 100
            if (i1 == 0) then
               i1 = conect(i,2)
               if (i1 == 0) GOTO 100
            end if
            a1 = atoms(i1)
            if (a1 /= iSilicon) GOTO 100
            ! if ( any(conect(i,2:) /= 0) ) GOTO 100
            if ( any(conect(i,3:) /= 0) ) GOTO 100
         case default
            stop 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ierr = i
   END SUBROUTINE


   SUBROUTINE check_proximity2_H(ifirst,ilast,atoms,conect,consistent,ierr)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: atoms(:),conect(:,:)
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k,ia,a1,a2,i1,i2
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         ia = atoms(i)
         select case(ia)
         case(iHydrogen)  ! Check H is bonded to 1 OH Oxygen
            k = conect(i,1)
            if (k == 0) GOTO 100
            if (atoms(k) /= iOxygenH) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
         case(iSilicon)  ! Check Si is bonded to 4 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100   ! but not Si
               if (atoms(k) == iHydrogen) GOTO 100  ! and not H
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            if (a1 /= iSilicon) GOTO 100
            i2 = conect(i,2)  !  must be bonded to  2 Silicons
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (a2 /= iSilicon) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
         case(iOxygenH)  ! OH Oxygen must be bonded to Si and H
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            i2 = conect(i,2)
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (i2 == 0) GOTO 100
            if (.NOT.((a1 == iHydrogen .and. a2 == iSilicon).or. &
                      (a2 == iHydrogen .and. a1 == iSilicon))) GOTO 100
            if ( any(conect(i,3:) /= 0) ) GOTO 100
         case default
            stop 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ierr = i
      RETURN
   END SUBROUTINE

END MODULE check_structure_mod


!>>INCLUDE 'readline_mod.f90'
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

!!>INCLUDE 'matvec3d_mod.f90'

MODULE matvec3d_mod
   USE precision_mod, only: wp
   implicit none
   public:: dotp_3d,len_3d,crossp_3d
   public:: transpose_3d,matvec_3d,matmul_3d,getinv3d
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
      real(wp):: transpose_3d(3,3)
      real(wp),intent(in):: m(3,3)
      integer:: i,j
      forall(i=1:3,j=1:3) transpose_3d(j,i) = m(i,j)
   END FUNCTION transpose_3d

   PURE FUNCTION matvec_3d(m,v)
      real(wp):: matvec_3d(3)
      real(wp),intent(in):: m(3,3)
      real(wp),intent(in):: v(3)
      matvec_3d(1) = m(1,1)*v(1) + m(1,2)*v(2) + m(1,3)*v(3)
      matvec_3d(2) = m(2,1)*v(1) + m(2,2)*v(2) + m(2,3)*v(3)
      matvec_3d(3) = m(3,1)*v(1) + m(3,2)*v(2) + m(3,3)*v(3)
   END FUNCTION matvec_3d

   PURE FUNCTION matmul_3d(m1,m2)
      real(wp),dimension(3,3):: matmul_3d
      real(wp),dimension(3,3),intent(in):: m1,m2
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
      real(wp),intent(in):: hmat(3,3)
      real(wp),intent(out):: hmati(3,3)
      real(wp):: det,odet
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


!>> INCLUDE 'pdb_read_mod.f90'
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

!>>INCLUDE 'bond_list_mod.f90'

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

END MODULE bond_list_mod

!>>INCLUDE 'crossbond_mod.f90'
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
      integer:: ifirst,neigh,ncell(125)
      real(wp):: r3(3),bond_lens(4),dr(3)
      logical:: consistent

      !  need to specify the changes in the box length directions for
      !  initialisation of the neighbour list.
      boxl(1:3) = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
      boxl2(1:3) = boxl(1:3)/2.0_wp

      call  DELETE_NLIST()

      print*, 'nix = ',nix
      print*, 'niy = ',niy
      print*, 'niz = ',niz

      call INIT_NLIST(-boxl2,boxl2,5.0_wp,natom)
      call NEW_NLIST(1,natom)

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


!!>INCLUDE 'command_line_mod.f90'

MODULE COMMAND_LINE_MOD
   implicit none
CONTAINS

   FUNCTION command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   END FUNCTION

   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,status)
      integer         , INTENT(IN)           :: NUMBER
      character(len=*), INTENT(OUT), OPTIONAL:: VALUE
      integer         , INTENT(OUT), OPTIONAL:: LENGTH
      integer         , INTENT(OUT), OPTIONAL:: status
      character(len=1000):: TMPVAL
      integer,EXTERNAL:: IARGC
      if (NUMBER < 0) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(status)) STATUS = 1
          RETURN
      else if (NUMBER > IARGC()) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(status)) STATUS = 2
          RETURN
      end if
      if (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
      if (PRESENT(LENGTH)) then
          if (PRESENT(VALUE)) then
              LENGTH = LEN_TRIM(VALUE)
          else
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          end if
      end if
      if (PRESENT(status)) STATUS = 0
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
      USE check_structure_mod

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
      integer:: nimage,natom_orig,neigh,ncell(125),iOx
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
         call read_conn_pdb(iu)
      CLOSE(iu)

      natom_max = 2*natom

      call prox_org


      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*,'pre-crossbond'
      print*, 'not consistent-wrong atoms bonded'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if

      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if

      boxl_0(1:3)=boxl(1:3)
      boxl2(1:3)=boxl(1:3)/2.0_wp
      natom_max = 2*natom
      allocate(del_list(2*natom))
      del_list=0
      call INIT_NLIST(-boxl2,boxl2,5.0_wp/angstrom,natom_max)
      call new_nlist(1,natom_max)
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
                  if(atom(j)==iSilicon) then
                     dr(1:3) = ps_rxyz(i,1:3) - rxyz(j,1:3)
                     call pbc(dr)
                     if ( dot_product(dr,dr) < cutoff2 )  then
                        dl_logical(j)=.TRUE.
                     end if
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

      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*,'pre-crossbond'
      print*, 'not consistent-wrong atoms bonded'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if

      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent-wrong atoms bonded - post crossbond'
      print'(/)'
      call write_atom_info(6,ifirst)
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

! Find the largest cluster of atoms
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print '(/a,i8)','n_cluster = ',n_cluster
      call cluster_count(natom,n_cluster,atomL,MaxCluster,ClusterCount)

      print '(/a,i8)','n_cluster = ',n_cluster
      print *,'MaxCluster = ',MaxCluster
      print *,'ClusterCount(MaxCluster) = ',ClusterCount(MaxCluster)
      print *,'sum(ClusterCount) = ',sum(ClusterCount)
      print *,'natom = ',natom

! Delete all atoms not in Largest cluster
      j = natom
      do
         if (atomL(j) /= MaxCluster) then
            atomL(j) = atomL(natom)
            call delete_atom(j)
         end if
         j = j - 1
         if (j == 0) exit
      end do

! convert singly bonded oxygens to OH
      do i = 1,natom
         if (atom(i) == iOxygen) then
            if (count(proximity(i,:) > 0) < 2) atom(i) = iOxygenH
         end if
      end do

! Delete all Si(OH)3 groups
      i = natom
      do
!print '(5i7)',i
!call write_atom_info(6,i)
         if (atom(i) == iSilicon) then
         if (noh(i) == 3) then
            call find_1st_atomtype(i,iOxygen,iOx)
            atom(iOx) = iOxygenH
            call delete_group0(i,iOx)
            i = natom + 1
         end if
         endif
         i = i - 1
         if (i <= 0) exit
      end do

!      ! remove floating atoms
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!
!     ! sort and delete atoms
!      call shell(floater,del_list)
!      do i = floater,1,-1
!         call delete_atom(del_list(i))
!      end do
!
!print*, 'floater = ', floater
!print*, 'natom pre-float del',natom
!
!      call re_init_nlist(boxl,5.0_wp)
!      call new_nlist
!
!      !! count number of silicons which are not 4 coordinated
!      !! and which have no OH group
!      ucoord1 = 0
!      ucoord2 = 0
!      ucoord3 = 0
!      coord_loop: do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!            do j=1,4
!               if (proximity(i,j)/=0) then
!                  if (atom(proximity(i,j))==iOxygenH) CYCLE coord_loop
!               end if
!            end do
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               ucoord1 = ucoord1 + 1
!            case(2)
!               ucoord2 = ucoord2 + 1
!            case(3)
!               ucoord3 = ucoord3 + 1
!            end select
!         end if
!         end if
!      end do coord_loop
!      print*, 'ucoord (Si_0_1) = ', ucoord1
!      print*, 'ucoord (Si_0_2) = ', ucoord2
!      print*, 'ucoord (Si_0_3) = ', ucoord3
!      print'(2/)'
!
!      call prox_org
!
!do
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !                                                      !
!      !              remove any 1 coord Si                   !
!      !              remove any 2 coord Si                   !
!      !           add OH to any 3 coord Si                   !
!      !                                                      !
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ndel = 0
!      del_list = 0
!      nat_temp = natom
!      at_loop: do i = 1,nat_temp
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               do j=1,ncmax(atom(i))
!                  if (proximity(i,j)==0) cycle
!                  if (atom(proximity(i,j))==iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!               ndel = ndel + 1
!               del_list(ndel) = i
!
!            case(2)
!               do j=1,ncmax(atom(i))
!                  if (proximity(i,j)==0) cycle
!                  if (atom(proximity(i,j))==iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!               ndel = ndel + 1
!               del_list(ndel) = i
!
!            case(3)
!               ! add on an OH to the empty spot
!               natom = natom + 1
!               ra(1:3) = rxyz(proximity(i,1),:)
!               rb(1:3) = rxyz(proximity(i,2),:)
!               rc(1:3) = rxyz(proximity(i,3),:)
!               r_ba = ra - rb
!               call pbc(r_ba)
!               r_ca = ra - rc
!               call pbc(r_ca)
!               mag_crossp_3d = sqrt(dot_product(crossp_3d(r_ba,r_ca),crossp_3d(r_ba,r_ca)))
!               temp = crossp_3d(r_ba,r_ca)/mag_crossp_3d
!               rxyz(natom,:) = rxyz(i,:) + temp*bondl_SiO
!               call pbc(rxyz(natom,:))
!               atom(natom) = iOxygenH
!               proximity(natom,1) = i
!               proximity(i,4) = natom
!            end select
!         end if
!         end if
!      end do at_loop
!
!      call shell(ndel,del_list)
!
!print*, 'ndel = ', ndel
!print*, 'natom pre-treat - del',natom
!
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-treat - del',natom
!print'(2/)'
!
!      call prox_org
!
!      ! remove floating atoms
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!
!     ! sort and delete atoms
!      call shell(floater,del_list)
!      do i = floater,1,-1
!         call delete_atom(del_list(i))
!      end do
!
!print*, 'floater = ', floater
!
!      CALL Init_HKNonLattice(natom)
!      CALL HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
!      del_list=0
!      ndel=0
!      do i=1,natom
!         if (atomL(i)/=MaxCluster) then
!            ndel=ndel+1
!            del_list(ndel)=i
!         end if
!      end do
!
!      call shell(ndel,del_list)
!      do i=ndel,1,-1
!         call delete_atom(del_list(i))
!      end do
!      deallocate(ClusterCount)
!
!
!      !! count number of silicons which are not 4 coordinated
!      !! and which have no OH group
!      ucoord1 = 0
!      ucoord2 = 0
!      ucoord3 = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!!            if (ANY(atom(proximity(i,:))==iOxygenH)) CYCLE
!            do j=1,4
!               if (proximity(i,j)/=0) then
!                  if (atom(proximity(i,j))==iOxygenH) CYCLE
!               end if
!            end do
!
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               ucoord1 = ucoord1 + 1
!            case(2)
!               ucoord2 = ucoord2 + 1
!            case(3)
!               ucoord3 = ucoord3 + 1
!            end select
!         end if
!         end if
!      end do
!print*, 'ucoord (Si_0_1) = ', ucoord1
!print*, 'ucoord (Si_0_2) = ', ucoord2
!print*, 'ucoord (Si_0_3) = ', ucoord3
!print'(2/)'
!
!      if (ucoord1==0 .AND. ucoord2==0 .AND. ucoord3==0) EXIT
!end do
!
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!
!print*, 'floater = ', floater
!print*, 'natom pre-float del',natom
!
!      do k=floater,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-float del',natom
!
!     call prox_org
!
!     ! go over and delete any clusters of floating atoms
!
!      CALL Init_HKNonLattice(natom)
!      CALL HKNonLattice(natom,proximity,n_cluster,atomL)
!      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
!      del_list=0
!      ndel=0
!      do i=natom,1-1
!         if (atomL(i)/=MaxCluster) then
!            ndel=ndel+1
!            del_list(ndel)=i
!         end if
!      end do
!
!      do i=1,ndel
!         call delete_atom(del_list(i))
!      end do
!
!do
!      !! delete any FULLY bonded si atoms which have an 3 OH's
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 0) then
!            select case(count(atom(proximity(i,:))==iOxygenH))
!            case(1)
!               ! do nothing
!            case(2)
!               ! do nothing
!            case(3)
!            ndel = ndel + 1
!            del_list(ndel) = i
!
!            do j = 1,ncmax(atom(i))
!               if (proximity(i,j)==0) stop 'screwed up'
!               if (atom(proximity(i,j)) == iOxygenH) then
!                  ndel = ndel + 1
!                  del_list(ndel) = proximity(i,j)
!               else if (atom(proximity(i,j)) == iOxygen) then
!                  atom(proximity(i,j)) = iOxygenH
!               end if
!            end do
!            case default
!               cycle
!            end select
!         end if
!         end if
!      end do
!
!      call shell(ndel,del_list)
!
!      do k=ndel,1,-1
!        call delete_atom(del_list(k))
!      end do
!
!      !! delete any TRIPLE bonded si atoms which have an 3 OH's
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 1) then
!            ndel = ndel + 1
!            del_list(ndel) = i
!            select case(count(atom(proximity(i,:))==iOxygenH))
!            case(1)
!               ! do nothing
!            case(2)
!               ! no nothing
!            case(3)
!               do j = 1,ncmax(atom(i))
!                  if (proximity(i,j)==0) CYCLE
!                  if (atom(proximity(i,j)) == iOxygenH) then
!                     ndel = ndel + 1
!                     del_list(ndel) = proximity(i,j)
!                  else if (atom(proximity(i,j)) == iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!            end select
!
!         end if
!         end if
!      end do
!
!      call shell(ndel,del_list)
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-triple - del',natom
!
!      !! convert all singly bonded oxygens to
!      do i = 1,natom
!         if (atom(i)/=iOxygen) CYCLE
!         if (count(proximity(i,:)==0) == 3) then
!            atom(i) = iOxygenH
!         end if
!      end do
!
!print*, 'post convert'
!
!      !! delete any double bonded si atoms which have an OH
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 2) then
!            do k=1,4
!               if (proximity(i,k)/=0) then
!                  if (atom(proximity(i,k))==iOxygenH) then
!                     ndel = ndel + 1
!                     del_list(ndel) = i
!                     do j = 1,ncmax(atom(i))
!                        if (proximity(i,j)==0) cycle
!                        if (atom(proximity(i,j)) == iOxygenH) then
!                           ndel = ndel + 1
!                           del_list(ndel) = proximity(i,j)
!                        else if (atom(proximity(i,j)) == iOxygen) then
!                           atom(proximity(i,j)) = iOxygenH
!                        end if
!                     end do
!                  end if
!               end if
!            end do
!         end if
!         end if
!      end do
!
!      call shell(ndel,del_list)
!
!print*, 'ndel = ', ndel
!print*, 'natom double-del',natom
!
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!      !! delete any singly bonded si atoms
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i)==iSilicon) then
!         if (count(proximity(i,:)==0) == 3) then
!            ndel = ndel + 1
!            del_list(ndel) = i
!            do j = 1,ncmax(atom(i))
!               if (proximity(i,j)==0) cycle
!               if (atom(proximity(i,j)) == iOxygen) then
!                  atom(proximity(i,j)) = iOxygenH
!               end if
!            end do
!         end if
!         end if
!      end do
!
!      !! analyse the bonds in the system
!      nSi_1_OH = 0
!      nSi_1_OH2 = 0
!      nSi_1_OH3 = 0
!      nSi_2_OH = 0
!      nSi_2_OH2 = 0
!      nSi_2_OH3 = 0
!      nSi_3_OH = 0
!      nSi_3_OH2 = 0
!      nSi_3_OH3 = 0
!      nSi_4_OH = 0
!      nSi_4_OH2 = 0
!      nSi_4_OH3 = 0
!      do i=1,natom
!         if (atom(i) == iSilicon) then
!            n_OH = 0
!            do j = 1,ncmax(atom(i))
!               if (atom(proximity(i,j)) == iOxygenH) then
!                  n_OH = n_OH + 1
!               end if
!            end do
!            coord_type = count(proximity(i,:)/=0)
!            select case(coord_type)
!            case(1)
!               select case(n_OH)
!               case(1)
!                  nSi_1_OH = nSi_1_OH + 1
!               case(2)
!                  nSi_1_OH2 = nSi_1_OH2 + 1
!               case(3)
!                  nSi_1_OH3 = nSi_1_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(2)
!               select case(n_OH)
!               case(1)
!                  nSi_2_OH = nSi_2_OH + 1
!               case(2)
!                  nSi_2_OH2 = nSi_2_OH2 + 1
!               case(3)
!                  nSi_2_OH3 = nSi_2_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(3)
!               select case(n_OH)
!               case(1)
!                  nSi_3_OH = nSi_3_OH + 1
!               case(2)
!                  nSi_3_OH2 = nSi_3_OH2 + 1
!               case(3)
!                  nSi_3_OH3 = nSi_3_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(4)
!               select case(n_OH)
!               case(1)
!                  nSi_4_OH = nSi_4_OH + 1
!               case(2)
!                  nSi_4_OH2 = nSi_4_OH2 + 1
!               case(3)
!                  nSi_4_OH3 = nSi_4_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            end select
!         else if (atom(i) == iOxygen) then
!            ! do nothing
!         else if (atom(i) == iOxygenH) then
!            ! do nothing
!         else
!            print*, 'error in atom types'
!            print*, 'atom(',i,') = ', atom_name(atom(i))
!            print'(2/)'
!            print*, 'program terminated prematurely'
!            print'(/)'
!            stop
!         end if
!      end do
!
!      !! bonding summary
!      print'(2/)'
!!      print*, 'nSi_0H   = ', nSi_OH
!!      print*, 'nSi_0H2  = ', nSi_OH2
!!      print*, 'nSi_0H3  = ', nSi_OH3
!      print*, 'Singly bonded silicons'
!      print*, 'nSi_1_OH  = ',  nSi_1_OH
!      print*, 'nSi_1_OH2 = ',  nSi_1_OH2
!      print*, 'nSi_1_OH3 = ',  nSi_1_OH3
!
!      print'(2/)'
!      print*, 'Double bonded silicons'
!      print*, 'nSi_2_OH  = ',  nSi_2_OH
!      print*, 'nSi_2_OH2 = ',  nSi_2_OH2
!      print*, 'nSi_2_OH3 = ',  nSi_2_OH3
!
!      print'(2/)'
!      print*, 'Triple bonded silicons'
!      print*, 'nSi_3_OH  = ',  nSi_3_OH
!      print*, 'nSi_3_OH2 = ',  nSi_3_OH2
!      print*, 'nSi_3_OH3 = ',  nSi_3_OH3
!
!      print'(2/)'
!      print*, 'Fully bonded silicons'
!      print*, 'nSi_4_OH  = ',  nSi_4_OH
!      print*, 'nSi_4_OH2 = ',  nSi_4_OH2
!      print*, 'nSi_4_OH3 = ',  nSi_4_OH3
!      print'(2/)'
!
!      ! if all ducks are in a row break out
!      if (nSi_1_OH==0 .AND.&
!          nSi_1_OH2==0 .AND.&
!          nSi_1_OH3==0 .AND.&
!          nSi_2_OH==0 .AND.&
!          nSi_2_OH2==0 .AND.&
!          nSi_2_OH3==0 .AND.&
!          nSi_3_OH==0 .AND.&
!          nSi_3_OH2==0 .AND.&
!          nSi_3_OH3==0 .AND.&
!          nSi_4_OH3==0) exit
!end do
!
!
!      !! convert all singly bonded oxygens to
!      do i = 1,natom
!         if (atom(i)/=iOxygen) CYCLE
!         if (count(proximity(i,:)==0) == 3) then
!            atom(i) = iOxygenH
!         end if
!      end do
!
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!
!   print*, 'floater = ', floater
!   print*, 'natom pre-float del',natom
!
!      do k=floater,1,-1
!         call delete_atom(del_list(k))
!      end do


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
      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
         print*, 'not consistent 2 - end'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if

      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
         print*, 'not consistent 1 - end'
      print'(/)'
      call write_atom_info(6,ifirst)
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
         write(iu,'(5i7)') i,(proximity(i,j),j = 1,ncmax(atom(i)))
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

   SUBROUTINE read_conn_pdb(iu)
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

   END SUBROUTINE read_conn_pdb


END PROGRAM etch_amorphous

