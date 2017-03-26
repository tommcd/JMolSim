!!>include 'precision_mod.f90'

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


!!>include 'math_const_mod.f90'

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


!!>include 'phys_const_mod.f90'

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
      real(wp),parameter:: ev         =  eplus                ! ev       j
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
!
END MODULE phys_const_mod


!!>include 'commandline_mod.f90'

MODULE COMMANDLINE_MOD
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

END MODULE COMMANDLINE_MOD

!!>include 'command_arg_mod.f90'

MODULE command_arg_mod
   USE precision_mod
   USE commandline_mod
   implicit none
   PRIVATE
   integer:: icur_arg = -1
   !
   ! Gets a command argument and assigns it to a variable.
   ! It checks the status of the argument and writes a record to stdout.
   !
   public:: next_command_argument
   interface next_command_argument
      module procedure next_command_argument_int, next_command_argument_real, &
                       next_command_argument_logical,next_command_argument_string
   end interface
!
CONTAINS

   SUBROUTINE next_command_argument_int(var,var_name,stat)
      integer,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_int

   SUBROUTINE next_command_argument_real(var,var_name,stat)
      real(wp),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_real

   SUBROUTINE next_command_argument_logical(var,var_name,stat)
      logical,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_logical

   SUBROUTINE next_command_argument_string(var,var_name,stat)
      character(*),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
      end if
      var = ''
      var(1:nc) = carg(1:nc)
      write (*,*) var_name,' = ',trim(var)
   END SUBROUTINE next_command_argument_string

END MODULE command_arg_mod


!!>include 'files_mod.f90'

MODULE files_mod
   public get_free_file_unit, myflush
CONTAINS

   SUBROUTINE get_free_file_unit(iu)
      implicit none
      character(len=*),parameter:: sub_name = 'get_free_file_unit'
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

!!>include 'rand_mod.f90'

MODULE rand_mod
! a random number generator (Numerical Recipies)
! with auxillary subroutines for storing and setting the
! state space
   USE precision_mod, only: wp,i4b
   implicit none
   private
   public:: rand,get_rand_state,set_rand_state,read_rand_state,write_rand_state,gasdev,randvec
   integer(i4b),private,save:: ix = -1, iy = -1
   integer(i4b),public,save:: irand_seed = -1
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

   FUNCTION gasdev()
      real(wp):: gasdev
      integer,SAVE:: iset = 0
      real(wp),SAVE:: gset
      real(wp):: fac,rsq,v1,v2
      if (iset == 0) then
111     v1 = 2.0_wp*rand() - 1.0_wp
        v2 = 2.0_wp*rand() - 1.0_wp
        rsq = v1*v1 + v2*v2
        if (rsq >= 1.0_wp .or. rsq == 0.0_wp) GOTO 111
        fac = sqrt(-2.0_wp*log(rsq)/rsq)
        gset = v1*fac
        gasdev = v2*fac
        iset = 1
      else
        gasdev = gset
        iset = 0
      end if
      RETURN
   END FUNCTION gasdev

   FUNCTION randvec(NDIM)
      integer,intent(in):: NDIM
      real(wp):: randvec(NDIM)
      integer:: i
      do i = 1,NDIM
         randvec(i) = rand()
      end do
   END FUNCTION randvec

END MODULE rand_mod

!!>include 'sort_mod.f90'

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


!!>include 'global_vars_mod.f90'


MODULE global_vars_mod
   USE precision_mod
   USE phys_const_mod
   real(wp),parameter:: Angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: ulength = 1.0e-10_wp*Angstrom
   real(wp),parameter:: uenergy = eV
   real(wp),parameter:: upressure = uenergy/(ulength**3)
   real(wp),parameter:: torr2Pa = 101325.0_wp/760.0_wp
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = k_boltz/eV   ! 8.6173423E-5_wp
   !real(wp),parameter:: qstar = sqrt(elm_coupling/(uEnergy**2))
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp),parameter:: eVTokJmol = eV*navo/1000.0_wp
   real(wp),parameter:: kJmolToeV = 1.0_wp/eVTokJmol
   real(wp),parameter:: CaltokJ = 0.004186800000000046_wp
   real(wp):: kBoltzT
   integer:: natom, natom_max,nfixed,nseed = 0
   integer:: n_GAS,n_atom_GAS,n_O2 = 0,n_atom_o2
   integer:: n_N2 = 0,n_atom_n2,n_CO2 = 0,n_atom_co2,n_h2o=0,n_atom_h2o=0
   integer:: nattached, nrelax

END MODULE global_vars_mod

!!>include 'comvel_mod.f90'

MODULE comvel_mod
   USE precision_mod
   USE global_vars_mod
   implicit none
   interface comvel
      module procedure  comvel_vec, comvel_arr
   end interface
CONTAINS

   SUBROUTINE comvel_vec(natom,TEMP,vxyz,mass)
      USE rand_mod
!   *******************************************************************
!   ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
!   **                                                               **
!   ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (unit) MASS.**
!   ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
!   ** MOLECULES, AND NON-LINEAR MOLECULES.  Vector version          **
!   **                                                               **
!   *******************************************************************
!
! INITIAL VELOCITY FOR MOVING PARTICLE
      integer,intent(in):: natom
      real(wp),intent(in):: mass(:,:),TEMP
      real(wp),intent(inout):: vxyz(:,:)
      real(wp):: RTEMP,SUMX,SUMY,SUMZ,sum_mass
      integer:: i
      RTEMP = sqrt(TEMP*K_ev)
      if (natom <= 0) RETURN
      do  i = 1, natom
         ! INITIAL VELOCITY FROM MAXWELLIAN DISTRIBUTION
         vxyz(i,1) = RTEMP*gasdev()
         vxyz(i,2) = RTEMP*gasdev()
         vxyz(i,3) = RTEMP*gasdev()
      end do
      vxyz = vxyz/sqrt(mass)
! REMOVE NET MOMENTUM
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      sum_mass = 0.0_wp
      do i = 1, natom
         SUMX = SUMX + mass(i,1)*vxyz(i,1)
         SUMY = SUMY + mass(i,2)*vxyz(i,2)
         SUMZ = SUMZ + mass(i,3)*vxyz(i,3)
         sum_mass = sum_mass + mass(i,1)
      end do
      SUMX = SUMX/sum_mass
      SUMY = SUMY/sum_mass
      SUMZ = SUMZ/sum_mass
      do i = 1, natom
         vxyz(i,1) = vxyz(i,1) - SUMX
         vxyz(i,2) = vxyz(i,2) - SUMY
         vxyz(i,3) = vxyz(i,3) - SUMZ
      end do
   END SUBROUTINE comvel_vec


   SUBROUTINE comvel_arr(nmol,nat,TEMP,vxyz,mass)
      USE rand_mod
!   *******************************************************************
!   ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
!   **                                                               **
!   ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (unit) MASS.**
!   ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
!   ** MOLECULES, AND NON-LINEAR MOLECULES. Array version            **
!   **                                                               **
!   *******************************************************************
!
! INITIAL VELOCITY FOR MOVING PARTICLE
      integer,intent(in):: nmol,nat
      real(wp),intent(in):: mass(:,:,:),TEMP
      real(wp),intent(inout):: vxyz(:,:,:)
      real(wp):: RTEMP,SUMX,SUMY,SUMZ,sum_mass
      integer:: i,j
      RTEMP = sqrt(TEMP*K_ev)
      if (nmol <= 0 .or. nat <= 0) RETURN
      do i = 1,nmol
      do j = 1,nat
         vxyz(i,j,1) = RTEMP*gasdev()
         vxyz(i,j,2) = RTEMP*gasdev()
         vxyz(i,j,3) = RTEMP*gasdev()
      end do
      end do
      vxyz = vxyz/sqrt(mass)
! REMOVE NET MOMENTUM
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      sum_mass = 0.0_wp
      do i = 1,nmol
      do j = 1,nat
         SUMX = SUMX + vxyz(i,j,1)*mass(i,j,1)
         SUMY = SUMY + vxyz(i,j,2)*mass(i,j,2)
         SUMZ = SUMZ + vxyz(i,j,3)*mass(i,j,3)
         sum_mass = sum_mass + mass(i,j,1)
      end do
      end do
      SUMX = SUMX/sum_mass
      SUMY = SUMY/sum_mass
      SUMZ = SUMZ/sum_mass
      do i = 1,nmol
      do j = 1,nat
         vxyz(i,j,1) = vxyz(i,j,1) - SUMX
         vxyz(i,j,2) = vxyz(i,j,2) - SUMY
         vxyz(i,j,3) = vxyz(i,j,3) - SUMZ
      end do
      end do
   END SUBROUTINE comvel_arr

END MODULE comvel_mod


!!>include 'grid_data_mod.f90'

MODULE grid_data_mod
   USE precision_mod
   USE global_vars_mod, only: Angstrom
   implicit none

   TYPE grid3D
      integer:: nbin(3)
      real(wp):: rlower(3)
      real(wp):: rupper(3)
      real(wp):: side(3)
      real(wp):: del(3)
      real(wp):: del2(3)
      real(wp):: deli(3)
      real(wp):: volbin
      real(wp),allocatable:: x(:)
      real(wp),allocatable:: y(:)
      real(wp),allocatable:: z(:)
   END TYPE grid3D

   INTERFACE new_grid3D_data
      module procedure new_grid3D_data3, new_grid3D_data4
   END INTERFACE new_grid3D_data

CONTAINS

   SUBROUTINE new_grid3D(rlower,rupper,dl,T)
      real(wp),intent(in):: rlower(3),rupper(3),dl
      type(grid3D),intent(out):: T
      integer:: i
      ! store grid data
      T%rlower = rlower
      T%rupper = rupper
      T%side = rupper-rlower
      T%nbin = nint(T%side/dl)
      T%del = T%side/T%nbin
      T%deli = 1.0_wp/T%del
      T%volbin = (T%del(1)*T%del(2)*T%del(3)*Angstrom**3)
      T%del2 = T%del*0.5_wp
      allocate(T%x(T%nbin(1)),T%y(T%nbin(2)),T%z(T%nbin(3)))
      forall(i=1:T%nbin(1)) T%x(i) = T%rlower(1) + (i-0.5_wp)*T%del(1)
      forall(i=1:T%nbin(2)) T%y(i) = T%rlower(2) + (i-0.5_wp)*T%del(2)
      forall(i=1:T%nbin(3)) T%z(i) = T%rlower(3) + (i-0.5_wp)*T%del(3)
   END SUBROUTINE new_grid3D


   SUBROUTINE new_grid3D_data3(T,dat)
      type(grid3D),intent(in):: T
      real(wp),allocatable,intent(out):: dat(:,:,:)
      ! Allocate array for grid data
      allocate( dat(T%nbin(1),T%nbin(2),T%nbin(3)) )
   END SUBROUTINE new_grid3D_data3


   SUBROUTINE new_grid3D_data4(T,nd4,dat)
      type(grid3D),intent(in):: T
      integer,intent(in):: nd4
      real(wp),allocatable,intent(out):: dat(:,:,:,:)
      ! Allocate array for grid data
      allocate( dat(T%nbin(1),T%nbin(2),T%nbin(3),nd4) )
   END SUBROUTINE new_grid3D_data4

   PURE FUNCTION R_GRID(T,ix,iy,iz)
      real(wp):: R_GRID(3)
      type(grid3D),intent(in):: T
      integer,intent(in):: ix,iy,iz
      R_GRID = (/ T%x(ix),T%y(iy),T%z(iz) /)
   END FUNCTION R_GRID

   PURE FUNCTION RU_GRID(T,ix,iy,iz)
      real(wp):: RU_GRID(3)
      type(grid3D),intent(in):: T
      integer,intent(in):: ix,iy,iz
      RU_GRID = (/ T%x(ix),T%y(iy),T%z(iz) /) + T%del2
   END FUNCTION RU_GRID

   PURE FUNCTION RL_GRID(T,ix,iy,iz)
      real(wp):: RL_GRID(3)
      type(grid3D),intent(in):: T
      integer,intent(in):: ix,iy,iz
      RL_GRID = (/ T%x(ix),T%y(iy),T%z(iz) /) - T%del2
   END FUNCTION RL_GRID

END MODULE grid_data_mod

!!>include 'list_mod.f90'

MODULE list_mod
   implicit none

!   integer,private,parameter:: nlist_max = 10000  ! static version
!   type mylist
!      integer:: n=0
!      integer:: ind(nlist_max)=0
!   end type mylist

   type mylist
      integer:: n=0
      integer,allocatable:: ind(:)
      integer:: nlist_max=0
   end type mylist

   interface assignment(=)
      module procedure vec2list
   end interface

   interface add_to_list
      module procedure add_int_to_list
      module procedure add_vec_to_list
      module procedure add_list_to_list
   end interface

   interface in_list
      module procedure in_vec
      module procedure in_lst
   end interface in_list

   interface operator(.in.)
      module procedure in_lst
   end interface

!   interface operator(.intersection.)
!      module procedure list_intersection
!   end interface

!   interface operator(.union.)
!      module procedure list_union
!   end interface

   PUBLIC :: ASSIGNMENT(=), OPERATOR(.in.)

CONTAINS

! comment ouf if mylist type is static
   SUBROUTINE new_list(L,N)
      type(mylist),intent(inout):: L
      integer,intent(in):: N
      if (allocated(L%ind)) deallocate(L%ind)
      allocate(L%ind(N))
      L%n = 0
      L%ind = 0
      L%nlist_max = N
   END SUBROUTINE new_list

! comment ouf if mylist type is static
   SUBROUTINE delete_list(L)
      type(mylist),intent(inout):: L
      if (allocated(L%ind)) deallocate(L%ind)
      L%n = 0
      L%nlist_max = 0
   END SUBROUTINE delete_list


   SUBROUTINE reset_list(L)
      type(mylist),intent(inout):: L
      L%n = 0
      L%ind = 0
   END SUBROUTINE reset_list


   SUBROUTINE add_int_to_list(lst,iat)
      type(mylist),intent(inout):: lst
      integer,intent(in):: iat
      lst%n = lst%n + 1
      !if (lst%n > nlist_max) STOP 'list is too large'  !static version
      if (lst%n > lst%nlist_max) STOP 'list is too large'
      lst%ind(lst%n) = iat
   END SUBROUTINE add_int_to_list

   SUBROUTINE add_vec_to_list(lst,vec)
      type(mylist),intent(inout):: lst
      integer,intent(in):: vec(:)
      integer:: nold
      nold = lst%n
      lst%n = lst%n + size(vec)
      !if (lst%n > nlist_max) STOP 'list is too large'  !static version
      if (lst%n > lst%nlist_max) STOP 'list is too large'
      lst%ind(nold+1:lst%n) = vec
   END SUBROUTINE add_vec_to_list

   SUBROUTINE add_list_to_list(lst,lst2)
      type(mylist),intent(inout):: lst
      type(mylist),intent(in):: lst2
      integer:: nold
      nold = lst%n
      lst%n = lst%n + lst2%n
      !if (lst%n > nlist_max) STOP 'list is too large'  !static version
      if (lst%n > lst%nlist_max) STOP 'list is too large'
      lst%ind(nold+1:lst%n) = lst2%ind(1:lst2%n)
   END SUBROUTINE add_list_to_list

   SUBROUTINE remove_from_list(lst,iat)
      type(mylist),intent(inout):: lst
      integer,intent(in):: iat
      integer:: i
      do i = 1,lst%n
         if (lst%ind(i) == iat) exit
      end do
      lst%ind(i) = lst%ind(lst%n)
      lst%ind(lst%n) = 0
      lst%n = lst%n - 1
   END SUBROUTINE remove_from_list


   SUBROUTINE remove_posn_list(lst,ipos)
      type(mylist),intent(inout):: lst
      integer,intent(in):: ipos
      lst%ind(ipos) = lst%ind(lst%n)
      lst%ind(lst%n) = 0
      lst%n = lst%n - 1
   END SUBROUTINE remove_posn_list


   SUBROUTINE rand_from_list(lst,iat)
      USE rand_mod
      integer,intent(out):: iat
      type(mylist),intent(inout):: lst
      integer:: j
      if (lst%n > 0) then
         j = int(rand()*lst%n) + 1
         iat = lst%ind(j)
      else
         iat = 0
      end if
   END SUBROUTINE rand_from_list


   SUBROUTINE vec2list(L,V)
      integer,intent(in):: V(:)
      type(mylist),intent(out):: L
      L%n = size(V)
      ! if (.not.allocated(L%ind)) allocate(L%ind(L%n))
      L%ind(1:L%n) = V
   END SUBROUTINE vec2list


   SUBROUTINE list_intersection(L1,L2,L)
      ! List L is the Intersection of L1 and L2
      ! L1 and L2 MUST BE SORTED & CONTAIN NO DUPLICATES
      type(mylist),intent(in):: L1,L2
      type(mylist),intent(inout):: L
      integer:: i,i1,i2,p
      i1 = 1
      i2 = 1
      i = 1
      L%ind = 0
      do
         if (i1 > L1%n .or. i2 > L2%n) exit
         p = L1%ind(i1)-L2%ind(i2)
         if (p < 0) then
            i1 = i1+1
         else if (p > 0) then
            i2 = i2+1
         else
            L%ind(i) = L1%ind(i1)
            i = i+1
            i1 = i1+1
            i2 = i2+1
         end if
      end do
      L%n = i-1
   END SUBROUTINE list_intersection


!   function intersection(L1,L2) result(L)
!      ! List L is the Intersection of L1 and L2
!      ! L1 and L2 MUST BE SORTED & CONTAIN NO DUPLICATES
!      type(mylist),intent(in):: L1,L2
!      type(mylist):: L
!      !allocate(L%ind(min(L1%n,L2%n)))
!      call list_intersect(L1,L2,L)
!   end function intersection


   SUBROUTINE list_union(L1,L2,L)
      ! List L is the Union of L1 and L2
      ! L1 and L2 MUST BE SORTED & CONTAIN NO DUPLICATES
      type(mylist),intent(in):: L1,L2
      type(mylist),intent(inout):: L
      integer:: i,i1,i2,p
      i1 = 1
      i2 = 1
      i = 1
      L%ind = 0
      do
         if (i1 > L1%n .or. i2 > L2%n) exit
         p = L1%ind(i1)-L2%ind(i2)
         if (p < 0) then  ! L1%ind(i1) < L2%ind(i2)
            L%ind(i) = L1%ind(i1)
            i = i + 1
            i1 = i1+1
         else if (p > 0) then ! L1 > L2
            L%ind(i) = L2%ind(i2)
            i = i + 1
            i2 = i2+1
         else  ! L2 == L1
            L%ind(i) = L1%ind(i1)
            i = i+1
            i1 = i1+1
            i2 = i2+1
         end if
      end do
      if (i1 > L1%n)then
         do
            if (i2 > L2%n) exit
            L%ind(i) = L2%ind(i2)
            i = i + 1
            i2 = i2+1
         end do
      else if(i2 > L2%n) then
         do
            if (i1 > L1%n) exit
            L%ind(i) = L1%ind(i1)
            i = i + 1
            i1 = i1+1
         end do
      end if
      L%n = i-1
   END SUBROUTINE list_union


   LOGICAL PURE FUNCTION in_vec(k,list,nl)
      integer,intent(in):: nl,list(:),k
      integer:: i
      in_vec = .FALSE.
      do i = 1,nl
         if (list(i) == k) then
            in_vec = .TRUE.
            RETURN
         end if
      end do
   END FUNCTION in_vec

   LOGICAL PURE FUNCTION in_lst(k,L)
      type(mylist),intent(in):: L
      integer,intent(in):: k
      integer:: i
      in_lst = .FALSE.
      do i = 1,L%n
         if (L%ind(i) == k) then
            in_lst = .TRUE.
            RETURN
         end if
      end do
   END FUNCTION in_lst


   SUBROUTINE rem_duplicates_list(nl,list)
      integer,intent(inout):: nl,list(:)
      integer:: i,j,a
! very crude ... just for testing
      i = 0
      do
         i = i + 1
         if (i > nl - 1) exit
         a = list(i)
         j = i
         do
            j = j + 1
            if (j > nl) exit
            if (a == list(j)) then
               list(j) = list(nl)
               nl = nl - 1
               j = j - 1
            end if
         end do
      end do
   END SUBROUTINE rem_duplicates_list

END MODULE list_mod

!!>include 'listset_mod.f90'

MODULE listset_mod
   implicit none

   integer,private,parameter:: nlist_max = 10000 ! Max number of elements
   integer,private,parameter::  max_card = 50000 ! Max element
   type listset
      integer:: n=0
      integer:: ind(nlist_max)=0
      integer:: element(max_card)=0
   end type listset
   type(listset),parameter:: null_listset = listset(0,0,0)

!   type listset
!      integer:: n=0
!      integer,allocatable:: ind(:)
!      integer,allocatable:: element(:)
!      integer:: nlist_max
!      integer:: max_card
!   end type listset

   interface assignment(=)
      module procedure vec2listset
   end interface

   interface add_to_list
      module procedure add_int_to_listset
      module procedure add_vec_to_listset
      module procedure add_list_to_listset
   end interface add_to_list

   interface remove_from_list
      module procedure remove_from_listset
   end interface remove_from_list

   interface remove_posn_list
      module procedure remove_posn_listset
   end interface remove_posn_list

   interface in_list
      module procedure in_listset
   end interface in_list

   interface operator(.in.)
      module procedure in_listset
   end interface

!   interface operator(.intersection.)
!      module procedure listset_intersection
!   end interface
!
!   interface operator(.union.)
!      module procedure listset_union
!   end interface

   PUBLIC:: ASSIGNMENT(=), OPERATOR(.in.)
!   PUBLIC:: operator(.intersection.), operator(.union.)

CONTAINS

!   SUBROUTINE new_listset(L,nlist_max,Max_Card)
!      type(listset),intent(inout):: L
!      integer,intent(in):: nlist_max,Max_Card
!      if (allocated(L%ind)) deallocate(L%ind,L%element)
!      allocate(L%ind(nlist_max),L%element(Max_Card))
!      L%n = 0
!      L%ind = 0
!      L%nlist_max = nlist_max
!      L%Max_Card = Max_Card
!   END SUBROUTINE new_listset
!
!   SUBROUTINE delete_listset(L)
!      type(listset),intent(inout):: L
!      if (allocated(L%ind)) deallocate(L%ind,L%element)
!      L%n = 0
!      L%nlist_max = 0
!      L%Max_Card = 0
!   END SUBROUTINE delete_listset


   SUBROUTINE add_int_to_listset(lst,iat)
      type(listset),intent(inout):: lst
      integer,intent(in):: iat
      if (lst%element(iat) > 0) RETURN ! already in list
      lst%n = lst%n + 1
      if (lst%n > nlist_max) STOP 'list is too large'
      if (iat > max_card) STOP 'element is too large'
      lst%ind(lst%n) = iat
      lst%element(iat) = lst%n
   END SUBROUTINE add_int_to_listset

   SUBROUTINE add_vec_to_listset(lst,iat)
      type(listset),intent(inout):: lst
      integer,intent(in):: iat(:)
      integer:: i
      do i = 1,size(iat)
         call add_int_to_listset(lst,iat(i))
      end do
   END SUBROUTINE add_vec_to_listset

   SUBROUTINE add_list_to_listset(lst,lst2)
      type(listset),intent(inout):: lst
      type(listset),intent(in):: lst2
      integer:: i
      do i = 1,lst2%n
         call add_int_to_listset(lst,lst2%ind(i))
      end do
   END SUBROUTINE add_list_to_listset


   SUBROUTINE remove_from_listset(lst,iat)
      type(listset),intent(inout):: lst
      integer,intent(in):: iat
      integer:: i
      i = lst%element(iat)
      lst%ind(i) = lst%ind(lst%n)
      lst%element(lst%ind(lst%n)) = i
      lst%ind(lst%n) = 0
      lst%n = lst%n - 1
      lst%element(iat) = 0
   END SUBROUTINE remove_from_listset


   SUBROUTINE remove_posn_listset(lst,ipos)
      type(listset),intent(inout):: lst
      integer,intent(in):: ipos
      integer:: j
      j = lst%ind(ipos)
      call remove_from_listset(lst,j)
   END SUBROUTINE remove_posn_listset


!   SUBROUTINE rand_from_listset(lst,iat)
!      USE rand_mod
!      integer,intent(out):: iat
!      type(listset),intent(inout):: lst
!      integer:: j
!      if (lst%n > 0) then
!         j = int(rand()*lst%n) + 1
!         iat = lst%ind(j)
!      else
!         iat = 0
!      end if
!   END SUBROUTINE rand_from_listset


   subroutine vec2listset(L,V)
      integer,intent(in):: V(:)
      type(listset),intent(out):: L
      integer:: i
      L = null_listset
      do i = 1,size(V)
         call add_int_to_listset(L,V(i))
      end do
   end subroutine vec2listset


   subroutine listset_intersection(L1,L2,L)
      ! List L is the Intersection of L1 and L2
      type(listset),intent(in):: L1,L2
      type(listset),intent(inout):: L
      integer:: i
      L = null_listset
      do i = 1, max(maxval(L1%ind(1:L1%n)),maxval(L2%ind(1:L2%n)))
         if ((L1%element(i) > 0) .and. (L2%element(i) > 0)) call add_int_to_listset(L,i)
      end do
   end subroutine listset_intersection

!   function listset_intersection(L1,L2) result(L)
!      ! List L is the Intersection of L1 and L2
!      type(listset),intent(in):: L1,L2
!      type(listset):: L
!      integer:: i
!      L = null_listset
!      do i = 1, max(maxval(L1%ind(1:L1%n)),maxval(L2%ind(1:L2%n)))
!         if ((L1%element(i) > 0) .and. (L2%element(i) > 0)) call add_int_to_listset(L,i)
!      end do
!   end function listset_intersection

   subroutine listset_union(L1,L2,L)
      ! List L is the Union of L1 and L2
      type(listset),intent(in):: L1,L2
      type(listset),intent(inout):: L
      integer:: i
      L = null_listset
      do i = 1, max(maxval(L1%ind(1:L1%n)),maxval(L2%ind(1:L2%n)))
         if ((L1%element(i) > 0) .or. (L2%element(i) > 0)) call add_int_to_listset(L,i)
      end do
   end subroutine listset_union

!   function listset_union(L1,L2) result(L)
!      ! List L is the Union of L1 and L2
!      type(listset),intent(in):: L1,L2
!      type(listset):: L
!      integer:: i
!      L = null_listset
!      do i = 1, max(maxval(L1%ind(1:L1%n)),maxval(L2%ind(1:L2%n)))
!         if ((L1%element(i) > 0) .or. (L2%element(i) > 0)) call add_int_to_listset(L,i)
!      end do
!   end function listset_union


   LOGICAL PURE FUNCTION in_listset(k,L)
      type(listset),intent(in):: L
      integer,intent(in):: k
      in_listset = (L%element(k) > 0)
   END FUNCTION in_listset

END MODULE listset_mod



!!>include 'HKNonLattice2_mod.f90'

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

!!>include 'unionfind_mod.f90'

module unionfind_mod
   implicit none
!
! Purpose: implements routines for a Union-Find alogithm for the purposes
!   of Hoshen-Kopelman cluster labeling
!
! Based on Tobin Fricke's implementation of the Hoshen-Kopelman algorithm for
! cluster labeling. Rewritten in FORTRAN by Jason Meiring
!
! Original notes:
!   Copyright (c) September 9, 2000, by Tobin Fricke <tobin@splorg.org>
!   Distributed under the terms of the GNU Public License.
!   Modified 2002-03-09 Tobin Fricke
!   Modified substantially 2004-04-21 by Tobin Fricke
!   http://splorg.org/people/tobin/projects/hoshenkopelman/hoshenkopelman.html
!
! Author: Jason Meiring
! Date: 2005/10/17
! Minor reformatting by T.C.McD 2006-10-5
!
! Implementation of Union-Find Algorithm:
!
! The 'labels' array has the meaning that labels(x) is an alias for the label
! x; by following this chain until x == labels(x), you can find the canonical
! name of an equivalence class.  The labels start at one; labels(0) is a
! special value indicating the highest label already used.

!  Module variables
   integer, allocatable, private :: labels(:)
   integer, private :: n_labels = 0           ! length of the labels array

contains

   function uf_find(x) result(y)
      ! Returns the canonical label for the equivalence class containing x
      integer :: y
      integer, intent(in) :: x
      integer :: w, z
      y = x
      w = x
      do
         if (labels(y) == y) exit
         y = labels(y)
      end do
      do
         if (labels(w) == w) exit
         z = labels(w)
         labels(w) = y
         w = z
      end do
   end function uf_find


   function uf_union(x,y) result (z)
      ! Joins two equivalence classes and returns the canonical label of the resulting class
      integer, intent(in) :: x, y
      integer :: z
      z = uf_find(y)
      labels(uf_find(x)) = z
   end function uf_union


   function uf_make_set() result (x)
      ! Creates a new equivalence class and returns its label
      integer :: x
      labels(0) = labels(0) + 1
      labels(labels(0)) = labels(0)
      x = labels(0)
   end function uf_make_set


   subroutine uf_initialize(max_labels)
      ! Sets up the data structures needed by the union-find implementation
      integer, intent(in) :: max_labels
      integer :: i, astat
      n_labels = max_labels
      allocate(labels(0:n_labels), stat=astat)
      if (astat /= 0) stop 'uf_initialize: unable to allocate memory for labels array'
      forall (i = 1:n_labels) labels(i) = i
      labels(0) = 0
   end subroutine uf_initialize


   subroutine uf_done()
      ! Frees the memory used by the union-find data structures
      n_labels = 0
      deallocate(labels)
   end subroutine uf_done

end module unionfind_mod

!!>include 'hk_lattice_uf_mod_old.f90'

module hk_lattice_uf_mod
!
! Hoshen-Kopelman for a 3D lattice
!
   implicit none
   public :: hk_lattice_uf

   type :: position
      integer :: x, y, z
   end type position

contains

   subroutine hk_lattice_uf(lattice,periodic,ncluster)
!
! Label clusters in the lattice by using the Hoshen-Kopelman cluster labeling algorithm
!
      use unionfind_mod
      integer,intent(inout):: lattice(:,:,:)  ! Holds cluster numbers
      logical,intent(in):: periodic(3)
      integer,intent(out):: ncluster
      integer,allocatable:: new_labels(:)    ! Holds label numbers
      integer:: i,astat,ix,iy,iz,Nx,Ny,Nz
      integer:: b,neighbors,maxlabel
      integer,parameter:: nmax=26   ! Maximum number of neighbors to check
      type(position):: nloc(nmax) ! Array of neighbor locations
!
      Nx = size(lattice,1)
      Ny = size(lattice,2)
      Nz = size(lattice,3)
!print*, 'Nx=',Nx,'Ny=',Ny,'Nz =',Nz 

      maxlabel = (Nx*Ny*Nz)/2
      allocate(new_labels(0:maxlabel), stat=astat)
      if (astat /= 0) stop 'hk_lattice_uf: unable to allocate memory for new_label array'

      ! Initialize the new_label array
      new_labels = 0

      ! Initialize the Union-Find algorithm
      call uf_initialize(maxlabel)

      ! Scan the lattice for connected cells
      do iz = 1,Nz
      do iy = 1,Ny
      do ix = 1,Nx

! print*,'ix =',ix,'iy =',iy,'iz= ',iz 

         ! Do not consider voids
         if (lattice(ix,iy,iz) == 0) cycle

         !
         ! Verify if each pore neighbor at directions -x, -y and -z is labeled
         !
         neighbors = 0

         if ((iz > 1).and.(ix > 1).and.(iy > 1)) then
            call add_neighbor(ix-1,iy-1,iz-1)
            call add_neighbor(ix,  iy-1,iz-1)
            call add_neighbor(ix+1,iy-1,iz-1)
!print*,'ix+1=',ix+1            
            call add_neighbor(ix-1,iy  ,iz-1)
            call add_neighbor(ix,  iy  ,iz-1)
            call add_neighbor(ix+1,iy  ,iz-1)
            call add_neighbor(ix-1,iy+1,iz-1)
            call add_neighbor(ix,  iy+1,iz-1)
            call add_neighbor(ix+1,iy+1,iz-1)
            call add_neighbor(ix-1,iy-1,iz)
            call add_neighbor(ix,  iy-1,iz)
            call add_neighbor(ix+1,iy-1,iz)
            call add_neighbor(ix-1,iy  ,iz)
         end if

!         if (periodic(3).and.(iz == Nz) then
!            call add_neighbor(ix,iy,1)
!         end if
!
!         if (periodic(2).and.(iy == Ny) then
!            
!            call add_neighbor(lattice(ix,1,iz)
!         end if
!
!         if (periodic(1)).and.(ix == Nx) then
!             call add_neighbor(lattice(1,iy-1,iz-1)
!             call add_neighbor(lattice(1,  iy,iz-1)
!             call add_neighbor(lattice(1,iy+1,iz-1)
!             call add_neighbor(lattice(1,iy-1,iz)
!             call add_neighbor(lattice(1,iy,iz)
!         end if

         ! Figure out what to do based on the number of bonded neighbors
         select case(neighbors)
         case(0)
            ! No neighbors -- give the cell a new cluster label
            lattice(ix,iy,iz) = uf_make_set()
         case(1)
            ! Part of an existing cluster. Give the current cell
            ! the same label as the cluster
            lattice(ix,iy,iz) = lattice(nloc(1)%x,nloc(1)%y,nloc(1)%z)
         case default
            ! More than one neighbor -- this site binds two or more clusters
            ! First, find the smallest label number and assign it to the
            ! current cell

            i = maxlabel+1
            do b = 1,neighbors
               if (lattice(nloc(b)%x,nloc(b)%y,nloc(b)%z) < i) then
                  i = lattice(nloc(b)%x,nloc(b)%y,nloc(b)%z)
               end if
            end do
            lattice(ix,iy,iz) = i

            ! Join the other clusters
            do b = 1,neighbors
               if (lattice(nloc(b)%x,nloc(b)%y,nloc(b)%z) /= i) then
                  lattice(ix,iy,iz) = uf_union(i,lattice(nloc(b)%x,nloc(b)%y,nloc(b)%z))
               end if
            end do
         end select

      end do
      end do
      end do

      ! Perform re-labeling
      do iz = 1,Nz
      do iy = 1,Ny
      do ix = 1,Nx
         if (lattice(ix,iy,iz) == 0) cycle
         i = uf_find(lattice(ix,iy,iz))
         if (new_labels(i) == 0) then
            new_labels(0) = new_labels(0) + 1;
            new_labels(i) = new_labels(0)
         end if
         lattice(ix,iy,iz) = new_labels(i)
      end do
      end do
      end do
      ncluster = new_labels(0)

      ! Free up memory
      deallocate(new_labels)
      call uf_done()
   contains
      subroutine add_neighbor(ix,iy,iz)
         integer,intent(in):: ix,iy,iz
         integer:: ixx,iyy,izz
         ixx=ix;iyy=iy;izz=iz
         !print*, 'ixx=',ixx
          if (periodic(1)) then
            if(ixx>nx)ixx=1
          end if
          if (periodic(2)) then
            if(iyy>ny)iyy=1
          end if
          if (periodic(3)) then
            if(izz>nz)izz=1
          end if
          !print*, 'ixx######=',ixx
          if(lattice(ixx,iyy,izz) /= 0) then
            neighbors = neighbors + 1
            nloc(neighbors) = position(ixx,iyy,izz)
          end if
      end subroutine add_neighbor
   end subroutine hk_lattice_uf



end module hk_lattice_uf_mod

!!>include 'percolation_lattice_mod.f90'

MODULE PERCOLATION_LATTICE_MOD
   USE precision_mod
   implicit none

   INTERFACE accessible_regions_z
      module procedure accessible_regions_z0
      module procedure accessible_regions_z1
   END INTERFACE accessible_regions_z

CONTAINS

   SUBROUTINE ACCESSIBLE_REGIONS_Z1(Lattice,nzL,nzU,Accessible,NumAcc)
      USE list_mod
      USE listset_mod
      !
      ! Given: a 3D array of labeled connected cells (clusters)
      ! where adjacent occupied cells have the same cluster label
      ! Return:
      ! (1) Accessible; the logical array of those cells between
      !     between layers z-layer = nzL to nzU
      !     where elements are marked true if they are
      !     'connected' to z-layer nzU and so have a cluster label
      !     that also ocurrs on z-layer nzU.
      ! (2) NumAcc; A list of the number of accessible cells in each
      !     z-layer from nzU down to the last connected layer
      !
      integer,intent(in):: Lattice(:,:,:),nzL,nzU
      logical,intent(out):: Accessible(:,:,:)
      type(mylist),intent(inout):: NumAcc
      type(listset):: L_top  ! Store the list of labels in layer nzU
      integer:: ix,iy,iz,lab,Nx,Ny,nacc,inc
      !
      Nx = size(Lattice,1)
      Ny = size(Lattice,2)

!print*,'nzl and nzu',nzL,nzU
!print*,'accessible_regions_z'
!print*,'Nx',Nx
!print*,'Ny',Ny
      !L_top = null_listset  ! Null set/ empty list
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,nzU)
         if (lab == 0) cycle ! Unoccupied 
         if ( .not.(lab.in.L_top) ) then
            call add_to_list(L_top,lab)
         end if
      end do
      end do

     Accessible = .FALSE.!original
      call reset_list(NumAcc)
      if (nzU > nzL) then
         inc = -1
      else
         inc = 1
      end if
      do iz = nzU,nzL,inc
 !print*,'iz=',iz     
         nacc = 0
         do iy = 1,Ny
         do ix = 1,Nx
            lab = Lattice(ix,iy,iz)
           if (lab == 0) cycle
           if ( lab.in.L_top ) then
              Accessible(ix,iy,iz) = .TRUE.
              nacc= nacc + 1
!print*,'ix = ',ix,'iy=',iy,'iz=',iz,'Acc',Accessible(ix,iy,iz) 
           end if
            
         end do
         end do
         if (nacc == 0) RETURN
         call add_to_list(NumAcc,nacc)
      end do
   END SUBROUTINE ACCESSIBLE_REGIONS_Z1


   SUBROUTINE ACCESSIBLE_REGIONS_Z0(Lattice,Accessible)
      USE list_mod
      USE listset_mod
      !
      ! Given: a 3D array of labeled connected cells (clusters)
      ! where adjacent occupied cells have the same cluster label
      ! Return:
      ! (1) Accessible; the logical array of those cells
      !     where elements are marked true if they are
      !     'connected' to the top z-layer (Nz) and so have a
      !     luster label that also ocurrs on z-layer Nz.
      !
      integer,intent(in):: Lattice(:,:,:)
      logical,intent(out):: Accessible(:,:,:)
      type(listset):: L_top  ! Store the list of labels in top layer
      integer:: ix,iy,iz,lab,Nx,Ny,Nz
      !
      Nx = size(Lattice,1)
      Ny = size(Lattice,2)
      Nz = size(Lattice,3)
 !print*,Nx,Ny,Nz
      !L_top = null_listset  ! Null set/ empty list
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,Nz)
         if (lab == 0) cycle ! Unoccupied
         if ( .not.(lab.in.L_top) ) then
            call add_to_list(L_top,lab)
         end if
      end do
      end do
      !
      do iz = 1,Nz
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,iz)
         if (lab == 0) cycle
         if ( lab.in.L_top ) then
            Accessible(ix,iy,iz) = .TRUE.
         else
            Accessible(ix,iy,iz) = .FALSE.
         end if
      end do
      end do
      end do
   END SUBROUTINE ACCESSIBLE_REGIONS_Z0


   FUNCTION percolation_check(lattice,dim_xyz)
      ! A crude check for percolation in dimension dim_xyz
      ! Returns the 1st spanning label found or -1 otherwise.
      integer,intent(in) :: lattice(:,:,:),dim_xyz
      integer:: percolation_check
      integer:: ix,iy,iz,ix2,iy2,iz2,Nx,Ny,Nz
!
      Nx = size(lattice,1)
      Ny = size(lattice,2)
      Nz = size(lattice,3)
      select case(dim_xyz)
      case(1)
         do iy = 1, Ny
         do iz = 1, Nz
         do iy2 = 1, Ny
         do iz2 = 1, Nz
            if (lattice(1,iy,iz)/= 0 .and. lattice(1,iy,iz) == lattice(Nx,iy2,iz2)) then
               percolation_check = lattice(1,iy,iz)
               return
            end if
         end do
         end do
         end do
         end do
      case(2)
         do ix = 1, Nx
         do iz = 1, Nz
         do ix2 = 1, Nx
         do iz2 = 1, Nz
            if (lattice(ix,1,iz)/= 0 .and. lattice(ix,1,iz) == lattice(ix2,Ny,iz2)) then
               percolation_check = lattice(ix,1,iz)
               return
            end if
         end do
         end do
         end do
         end do
      case(3)
         do ix = 1, Nx
         do iy = 1, Ny
         do ix2 = 1, Nx
         do iy2 = 1, Ny
            if (lattice(ix,iy,1)/= 0 .and. lattice(ix,iy,1) == lattice(ix2,iy2,Nz)) then
               percolation_check = lattice(ix,iy,1)
               return
            end if
         end do
         end do
         end do
         end do
      end select
      percolation_check = -1
   END FUNCTION percolation_check


   subroutine percolation_data_print(iu,array,perc_label)
      integer,intent(in) :: iu,array(:,:,:),perc_label
      integer:: ix,iy,iz
      if (perc_label == (-1)) then
         write(iu,*) "-1 There was no percolation"
      else
         do iz = 1,size(array,3)
         do iy = 1,size(array,2)
         do ix = 1,size(array,1)
            if (array(ix,iy,iz) == perc_label) write(iu,"('c',i0,3i6)") array(ix,iy,iz),ix,iy,iz
         end do
         end do
         end do
      end if
   end subroutine percolation_data_print


   subroutine print_mat(iu,array)
      integer,intent(in) :: iu,array(:,:,:)
      integer:: ix,iy,iz
      do iz = 1,size(array,3)
      write(iu,*) 'iz = ',iz
      do iy = 1,size(array,2)
      do ix = 1,size(array,1)
         if (array(ix,iy,iz) /= 0) then
            write(iu,"((1x,i5))",advance='no') array(ix,iy,iz)
         else
            write(iu,"(a6)",advance='no') '  __  '
         end if
      end do
      write(iu,*)
      end do
      write(iu,*)
      end do
   end subroutine print_mat


   subroutine init_lattice(threshold,array,lattice)
      use precision_mod
      real(wp),intent(in):: threshold
      real(wp),intent(in):: array(:,:,:)
      integer,intent(out):: lattice(:,:,:)
      integer:: ix,iy,iz
      do iz = 1,size(array,3)
      do iy = 1,size(array,2)
      do ix = 1,size(array,1)
         ! Put 1's in places < threshold
         if (array(ix,iy,iz) <= threshold) then
            lattice(ix,iy,iz) = 1
         else
            lattice(ix,iy,iz) = 0
         end if
      end do
      end do
      end do
   end subroutine init_lattice


   subroutine check_labelling(lattice)
     ! This procedure checks to see that any occupied neighbors of an occupied site
     ! have the same label.
      integer,intent(in):: lattice(:,:,:)
      integer,parameter:: nmax=26
      integer:: nclus(nmax),neigh
      integer:: i,ic,ix,iy,iz,Nx,Ny,Nz
      Nx = size(lattice,1)
      Ny = size(lattice,2)
      Nz = size(lattice,3)
      do iz = 1,Nz
      do iy = 1,Ny
      do ix = 1,Nx
         ic = lattice(ix,iy,iz)
         neigh = 0
         nclus = 0
         if (ic == 0) cycle
         ! list the neighbor labels
         if (ix == 1) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix-1,iy,iz))
         end if
         if (ix == Nx) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix+1,iy,iz))
         end if
         if (iy == 1) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix,iy-1,iz))
         end if
         if (iy == Ny) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix,iy+1,iz))
         end if
         if (iz == 1) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix,iy,iz-1))
         end if
         if (iz == Nz) then
            call add_neigh(0)
         else
            call add_neigh(lattice(ix,iy,iz+1))
         end if
         do i = 1, neigh
            if (nclus(i) == 0) cycle
            if (nclus(i) /= ic) then
               stop 'check_labelling: error'
            end if
         end do
      end do
      end do
      end do
   contains
      subroutine add_neigh(ii)
         integer,intent(in):: ii
         neigh = neigh + 1
         nclus(neigh) = ii
      end subroutine add_neigh
   end subroutine check_labelling


   subroutine statistics_lattice(lattice,ncluster, meanCL, stdDevCL, minCL, maxCL)
!
! Determine the number, average size, and std. deviation size of
! clusters in the lattice
!
      integer, intent(in) :: lattice(:,:,:)
      integer, intent(out) :: ncluster  ! Number of individual clusters
      real, intent(out) :: meanCL       ! mean clustersize
      real, intent(out) :: stdDevCL     ! std. deviation of clusters
      integer,intent(out) :: minCL      ! smallest cluster
      integer, intent(out) :: maxCL     ! largest cluster
!
      integer, allocatable :: sizes(:)        ! Holds clusters sizes
      integer :: i,astat,ix,iy,iz,Nx,Ny,Nz,maxlabel

      ! Allocate and initialize sizes array
      Nx = size(lattice,1)
      Ny = size(lattice,2)
      Nz = size(lattice,3)
      maxlabel = (Nx*Ny*Nz)/2
      allocate(sizes(1:maxlabel), stat=astat)
      if (astat /= 0) stop 'statistics_lattice: unable to allocate memory for sizes array'
      do i = 1,maxlabel
         sizes(i) = 0
      end do

      ! Compute the sizes (lengths) of the clusters based on their cluster number
      do iz = 1,Nz
      do iy = 1,Ny
      do ix = 1,Nx
         i = lattice(ix,iy,iz)
         if (i /= 0) then
            sizes(i) = sizes(i) + 1
         end if
      end do
      end do
      end do

      ! Compute statistics on clusters of size > 1
      minCL = maxlabel+1
      maxCL = 0
      meanCL = 0.0
      ncluster = 0

      do i = 1,maxlabel
         if (sizes(i) == 0) exit
         if (sizes(i) > 1) then
            ncluster = ncluster + 1
            if (sizes(i) < minCL) minCL = sizes(i)
            if (sizes(i) > maxCL) maxCL = sizes(i)
            meanCL = meanCL + sizes(i)
         end if
      end do

      if (ncluster > 0) then
         ! Compute the average cluster length
         meanCL = meanCL/real(ncluster)
         ! Compute the standard deviation of cluster length
         stdDevCL = 0.0
         do i = 1,maxlabel
            if (sizes(i) == 0) exit
            if (sizes(i) > 1) then
               stdDevCL = stdDevCL + (real(sizes(i)) - meanCL)**2
            end if
         end do
         stdDevCL = stdDevCL / real(ncluster)
         stdDevCL = sqrt(stdDevCL)
      else
         meanCL = 0.0
         stdDevCL = 0.0
         minCL = 0
         maxCL = 0
      endif

      ! Free up memory
      deallocate(sizes)
   end subroutine statistics_lattice


   SUBROUTINE List_Lattice_Span_Z(Lattice,L_span)
      USE listset_mod
      !
      ! Given: a 3D array of labeled connected cells (clusters)
      ! where adjacent occupied cells have the same cluster label
      ! Return:
      !    L_span; A list of the cluster numbers which span
      !            in the z-dimension
      !
      integer,intent(in):: Lattice(:,:,:)
      type(listset),intent(inout):: L_span
      type(listset):: L_top, L_bot
      integer:: ix,iy,lab,Nx,Ny,Nz
      !
      Nx = size(Lattice,1)
      Ny = size(Lattice,2)
      Nz = size(Lattice,3)
      ! List clusters in top layer
      !L_top = null_listset  ! Null set/ empty list
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,Nz)
         if (lab == 0) cycle ! Unoccupied
         if ( .not.(lab.in.L_top) ) then
            call add_to_list(L_top,lab)
         end if
      end do
      end do

      ! List clusters in bottom layer
      !L_bot = null_listset
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,1)
         if (lab == 0) cycle ! Unoccupied
         if ( .not.(lab.in.L_bot) ) then
            call add_to_list(L_bot,lab)
         end if
      end do
      end do

!      !L_span = null_listset
!      call listset_intersection(L_bot,L_top,L_span)
!
!print '(a,i6)','L_bot%n = ',L_bot%n
!print '(20i6)',L_bot%ind(1:L_bot%n)
!print '(a,i6)','L_top%n = ',L_top%n
!print '(20i6)',L_top%ind(1:L_top%n)
!print '(a,i6)','L_span%n = ',L_span%n
!print '(20i6)',L_span%ind(1:L_span%n)

      !L_span = null_listset
      do iy = 1,Ny
      do ix = 1,Nx
         lab = Lattice(ix,iy,1)
         if (lab == 0) cycle
         if (lab.in.L_top) then
            if (.not.(lab.in.L_span)) then
               call add_to_list(L_span,lab)
            end if
         end if
      end do
      end do
!print '(a,i6)','L_span%n = ',L_span%n
!print '(20i6)',L_span%ind(1:L_span%n)
   END SUBROUTINE List_Lattice_Span_Z

END MODULE PERCOLATION_LATTICE_MOD

!!>include 'Dreiding_parameters_mod.f90'


MODULE Dreiding_parameters_mod
!
! Parameters for the simplified Dreiding potential.
! Source: William A. Goddard, et al., J.Phys.chem, 94, 8897-8909(1990)
!
   USE precision_mod, only: wp
   USE global_vars_mod, only: angstrom, K_ev, qstar, pi
   USE phys_const_mod, only : conv4
   implicit none
! Geometric valence parameter (bond radius) from Dreiding
  real(wp),parameter:: rHc = 0.330_wp/Angstrom
  real(wp),parameter:: rOc = 0.560_wp/Angstrom
  real(wp),parameter:: rC2 = 0.670_wp/Angstrom
  real(wp),parameter:: rC3 = 0.770_wp/Angstrom
  real(wp),parameter:: rSi = 0.937_wp/Angstrom
  real(wp),parameter,private:: del = 0.010_wp/Angstrom

! Bond parameters
   real(wp),parameter:: K_bn_dreid = 30.352_wp*(Angstrom*Angstrom)  ! ((700 kcal/mol)/A2)eV A^-2
   real(wp),parameter:: ASiOc = (rSi + rOc - del)  ! 1.487_wp Angstrom
   real(wp),parameter:: AC2C3 = (rC2 + rC3 - del)  ! 1.430_wp Angstrom
   real(wp),parameter:: AOcC2 = (rOc + rC2 - del)  ! 1.320_wp Angstrom
   real(wp),parameter:: AC2Hc = (rC2 + rHc - del)  ! 1.090_wp Angstrom
   real(wp),parameter:: AC3Hc = (rC3 + rHc - del)  ! 1.090_wp Angstrom
   real(wp),parameter:: AOH   = (rOc + rHc - del)  ! 0.880_wp Angstrom
!
! Note:
!                      U_bond_Dreiding = (kd/2)*(d - d0)^2
! but we use this form  U_bond_Keating =  (k/2)*(d^2 - d0^2)^2
!
! =>  k ~ kd/(4*d0^2)
!
   real(wp),parameter:: KSiOc = (K_bn_dreid)/(4.0_wp*ASiOc**2)
   real(wp),parameter:: KC2C3 = (K_bn_dreid)/(4.0_wp*AC2C3**2)
   real(wp),parameter:: KOcC2 = (K_bn_dreid)/(4.0_wp*AOcC2**2)
   real(wp),parameter:: KC2Hc = (K_bn_dreid)/(4.0_wp*AC2Hc**2)
   real(wp),parameter:: KC3Hc = (K_bn_dreid)/(4.0_wp*AC3Hc**2)
   real(wp),parameter:: KOH   = (K_bn_dreid)/(4.0_wp*AOH**2)

!Angle parameters
   real(wp),parameter:: K_ang_dreid =  4.336_wp  ! ((100 kcal/mol)/rad2)
   real(wp),parameter:: cos_HCH   = -1.0_wp/3.0_wp
   real(wp),parameter:: sin_HCH   = 0.942809_wp
   real(wp),parameter:: cos_SiOC  = -0.5_wp
   real(wp),parameter:: sin_SiOC  = 0.866025_wp
   real(wp),parameter:: cos_OC2C3 = -1.0_wp/3.0_wp
   real(wp),parameter:: sin_OC2C3 = 0.942809_wp
!   real(wp),parameter:: angle_SiOH = (109.47122062658917679_wp/180.0_wp)*pi
   real(wp),parameter:: cos_SiOH  = -1.0_wp/3.0_wp
   real(wp),parameter:: sin_SiOH  =  0.942809_wp
!
! Note:
!                      U_angle_Dreiding = (kd/2)*(th - th0)^2
! but we use this form  U_angle_Keating =  (k/2)*(cos(th) - cos(th0))^2
!
! =>  k ~ kd/(sin(th0)^2)
!
   real(wp),parameter:: KOC2C3 =  (K_ang_dreid)/(sin_OC2C3**2)
   real(wp),parameter:: KSiOcC2 = (K_ang_dreid)/(sin_SiOC**2)
   real(wp),parameter:: KOcC2Hc = (K_ang_dreid)/(sin_OC2C3**2)
   real(wp),parameter:: KC2C3Hc = (K_ang_dreid)/(sin_OC2C3**2)
   real(wp),parameter:: KHcC3Hc = (K_ang_dreid)/(sin_OC2C3**2)
   real(wp),parameter:: KSiOH   = (K_ang_dreid)/(sin_SiOH**2)
   real(wp),parameter:: KHCH    = (K_ang_dreid)/(sin_HCH**2)
END MODULE Dreiding_parameters_mod

!!>include 'seaton_mod.f90'

MODULE seaton_mod
   USE precision_mod
   USE global_vars_mod, only: Angstrom, K_ev, qstar, pi
   implicit none
! atom types
   real(wp),parameter:: eps_O_sil    = 185.0_wp*K_ev,    sig_O_sil = 2.708_wp/Angstrom,  q_O_sil = -0.64025_wp*qstar
   real(wp),parameter:: eps_Si_sil   =   0.0_wp*K_ev,   sig_Si_sil = 0.0_wp/Angstrom,   q_Si_sil = -2.0_wp*q_O_sil
   real(wp),parameter:: eps_O_OH_sil = 185.0_wp*K_ev, sig_O_OH_sil = 3.0_wp/Angstrom, q_O_OH_sil = -0.533_wp*qstar
   real(wp),parameter:: eps_H_OH_sil =   0.0_wp*K_ev, sig_H_OH_sil = 0.0_wp/Angstrom, q_H_OH_sil =  0.206_wp*qstar
   real(wp),parameter:: eps_H_H2O   =  0.00_wp*K_ev,  sig_H_H2O = 0.00_wp/Angstrom,   q_H_H2O =  0.417_wp*qstar
   real(wp),parameter:: eps_O_H2O   = 76.58_wp*K_ev,  sig_O_H2O = 3.1506_wp/Angstrom, q_O_H2O = -2.0_wp*q_H_H2O
   real(wp),parameter:: eps_O_CO2   = 82.997_wp*K_ev,sig_O_CO2 = 3.064_wp/Angstrom, q_O_CO2 = -0.33225_wp*qstar
   real(wp),parameter:: eps_C_CO2   = 29.999_wp*K_ev,sig_C_CO2 = 2.785_wp/Angstrom, q_C_CO2 = -2.0_wp*q_O_CO2
   real(wp),parameter:: eps_n_N2    = 34.897_wp*K_ev, sig_n_N2 = 3.3211_wp/Angstrom, q_n_N2 = -0.5475_wp*qstar
   real(wp),parameter:: eps_O_O2    = 43.183_wp*K_ev, sig_O_O2 = 3.1062_wp/Angstrom, q_O_O2 = -0.3577_wp*qstar
   real(wp),parameter:: q_OH_sil =  q_O_OH_sil + q_H_OH_sil
!  Bond and Angle parameters
   real(wp),parameter:: bondl_SiO = 1.600_wp/Angstrom
   real(wp),parameter:: bondl_OH  = 0.945_wp/Angstrom
   real(wp),parameter:: bondl_O2  = 0.9699_wp/Angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/Angstrom
   real(wp),parameter:: bondl_CO2 = 1.1610_wp/Angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: theta_CO2 = pi
   real(wp),parameter:: K_CO2 = 5.331_wp   ! eV/Rad^2   k/2 (th - th0)^2
   real(wp),parameter:: KOO  = 3.66142_wp*(Angstrom*Angstrom)   ! eV A^-2  !  k/2 (d - d0)^2
   real(wp),parameter:: KNN  = 3.66142_wp*(Angstrom*Angstrom)   ! eV A^-2  !  k/2 (d - d0)^2
   real(wp),parameter:: KCO2 = 3.66142_wp*(Angstrom*Angstrom)   ! eV A^-2  !  k/2 (d - d0)^2
   real(wp),parameter:: AOO = bondl_O2
   real(wp),parameter:: ANN = bondl_N2
   real(wp),parameter:: ACO2 = bondl_CO2
END MODULE seaton_mod

!!>include 'charmm_mod.f90'

MODULE charmm_mod
   USE precision_mod
   USE global_vars_mod, only: Angstrom, K_ev, qstar, pi
   USE phys_const_mod, only : conv4
   implicit none
! LJ parameters & charges
   ! LJ(par_all_22_prot_inp), charge(top_all_22_prot_inp)
   real(wp),parameter:: eps_C2_sil = 0.055_wp*conv4*K_ev, sig_C2_sil = 3.0_wp/Angstrom, q_C2_sil = -0.18_wp*qstar
   real(wp),parameter:: eps_C3_sil = 0.080_wp*conv4*K_ev, sig_C3_sil = 3.0_wp/Angstrom, q_C3_sil = -0.27_wp*qstar
   real(wp),parameter::   eps_H_HC = 0.028_wp*conv4*K_ev,   sig_H_HC = 1.0_wp/Angstrom,   q_H_HC =  0.09_wp*qstar
   real(wp),parameter:: eps_Ch2_sil = 0.2150_wp*conv4*K_ev, sig_Ch2_sil = 3.6_wp/Angstrom, q_Ch2_sil = 0.0_wp*qstar !drieding
   real(wp),parameter:: eps_Ch3_sil = 0.3050_wp*conv4*K_ev, sig_Ch3_sil = 3.7_wp/Angstrom, q_Ch3_sil = 0.0_wp*qstar !drieding

   ! LJ(par_silicate_inp), charge(top_silicate_inp)
   real(wp),parameter:: eps_O_OC_sil = 185.0_wp*K_ev, sig_O_OC_sil = 2.708_wp/Angstrom, q_O_OC_sil = -0.64025_wp*qstar
   real(wp),parameter::   eps_Ox_sil = 0.0_wp*K_ev,     sig_Ox_sil = 0.0_wp/Angstrom
!  Bond and Angle parameters
   real(wp),parameter:: bondl_SiOc = 1.619_wp/Angstrom  !Teos(J.Phy.Chem.Vol.99,No.7,1995,2166 - 76)
   real(wp),parameter:: bondl_OC   = 1.407_wp/Angstrom  !Teos(J.Phy.Chem.Vol.99,No.7,1995,2166 - 76)
   real(wp),parameter:: bondl_CC   = 1.516_wp/Angstrom  !Teos(J.Phy.Chem.Vol.99,No.7,1995,2166 - 76)
   real(wp),parameter:: bondl_CH   = 1.083_wp/Angstrom  !Teos(J.Phy.Chem.Vol.99,No.7,1995,2166 - 76)
   real(wp),parameter:: angle_SiOC = (120.13_wp/180.0_wp)*pi !Methoxysilane(Acta.Crys. C44,1988,1-4)
END MODULE charmm_mod

!!>include 'Keating_parameters_mod_vonAlfthan.f90'

MODULE Keating_parameters_mod
!
! Parameters for the simplified Keating potential.
! Source: S. von alfthan et al., Phys. Rev. B 68, 073203, (2003)
!
   USE precision_mod, only: wp
   USE global_vars_mod, only: Angstrom
   implicit none
   real(wp),parameter:: KSiSi =  9.08_wp*(Angstrom**2)  ! eV A^-2
   real(wp),parameter:: KSiO  = 27.00_wp*(Angstrom**2)
   real(wp),parameter:: ASiSi = 2.35_wp/Angstrom
   real(wp),parameter:: ASiO =  1.61_wp/Angstrom
   real(wp),parameter:: KSiSiSi = 3.58_wp   ! eV
   real(wp),parameter:: KSiSiO  = 3.93263_wp
   real(wp),parameter:: KOSiO   = 4.32_wp
   real(wp),parameter:: KSiOSi  = 2.0_wp
   real(wp),parameter:: cos_OSiO  = -1.0_wp/3.0_wp
   real(wp),parameter:: cos_SiOSi = -0.809_wp
END MODULE Keating_parameters_mod



!!>include 'atom_types_mod.f90'
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
   integer,parameter:: iCarbon2 = 11
   integer,parameter:: iCarbon3 = 12
   integer,parameter:: iHydrogenC = 13
   integer,parameter:: iOxygenC = 14
   integer,parameter:: iOxygenF = 15
   integer,parameter:: iCarbonH2 = 16
   integer,parameter:: iCarbonH3 = 17
   integer,parameter:: ntyp = 17
   integer:: n_atomtype(ntyp)
   character(4),parameter:: atom_name(0:ntyp) = (/ ' __ ', &
      'Si  ','O   ','OH  ','H   ','Ow  ', &
      'Hw  ','O_CO','C_CO','N_N2','O_O2', &
      'C2  ','C3  ','Hc  ','Oc  ','Ox  ','Ch2 ','Ch3 ' /)
   integer,parameter:: ncmax(0:ntyp) = (/ 0, &
                       4, 2, 2, 1, 2, &
                       1, 2, 4, 3, 2, &
                       4, 4, 1, 2, 2, 2, 1 /)
   real(wp),parameter:: mSi = 28.06_wp
   real(wp),parameter:: mOx = 16.00_wp
   real(wp),parameter:: mOH = 17.00_wp
   real(wp),parameter:: mH =  1.00_wp
   real(wp),parameter:: mN = 14.00_wp
   real(wp),parameter:: mC = 12.00_wp
   real(wp),parameter:: mCh2 = 14.00_wp
   real(wp),parameter:: mCh3 = 15.00_wp
   real(wp),parameter:: amass(0:ntyp) = (/ 0.0_wp, &
                        mSi, mOx, mOH, mH, mOx, &
                         mH, mOx,  mC, mN, mOx, &
                         mC, mC,   mH, mOx,mOx,mCh2,mCh3 /)
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

!!>include 'teos_atomlst_mod.f90'

MODULE teos_atomlst_mod
   USE atom_types_mod
   implicit none
!
   integer,parameter:: L_OH     = 1
   integer,parameter:: L_OSiOH3 = 2
   integer,parameter:: L_OSiO   = 3
   integer,parameter:: L_SiOC3  = 4
   integer,parameter:: L_SiOHOC2 = 5
   integer,parameter:: L_SiOH2OC = 6
   integer,parameter:: L_SiOHOC2_SiO = 7
   integer,parameter:: L_SiOH2OC_SiO = 8
   integer,parameter:: L_SiOx3OC1 = 9
   integer,parameter:: L_SiOx2OC2 = 10
   integer,parameter:: L_SiOx2OCOH = 11 
   integer,parameter:: L_OC = 12

!
   integer,parameter:: iOSiOH3 = ntyp+1
   integer,parameter:: iOSiO   = ntyp+1
   integer,parameter:: iOHB    = ntyp+3
   integer,parameter:: iSiOC3 = ntyp + 4
   integer,parameter:: iSiOHOC2 = ntyp + 5
   integer,parameter:: iSiOH2OC = ntyp + 6
   integer,parameter:: iSiOHOC2_SiO = ntyp + 7
   integer,parameter:: iSiOH2OC_SiO = ntyp + 8
   integer,parameter:: iSiOx3OC1 = ntyp + 9
   integer,parameter:: iSiOx2OC2 = ntyp + 10
   integer,parameter:: iSiOx3OH1 = ntyp + 11
   integer,parameter:: iSiOx2OCOH = ntyp + 12
   integer,parameter:: iSiOx4 = ntyp + 13
   integer,parameter:: iSiOx2OH2 = ntyp + 14
   integer,parameter:: iSiOcOH3 = ntyp + 15
   integer,parameter:: natom_type_rxn = ntyp+15

END MODULE teos_atomlst_mod

!!>include 'constants_mod.f90'
MODULE constants_mod
   USE atom_types_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: Angstrom,pi
   USE keating_parameters_mod
   USE seaton_mod, only: koo, knn, kco2, aoo, ann, aco2, K_CO2, THETA_CO2, bondl_SiO, &
                         sig_O_OH_sil, eps_O_OH_sil, q_O_OH_sil,q_OH_sil, sig_O_sil, eps_O_sil,q_O_sil, &
                         sig_Si_sil, eps_Si_sil,q_Si_sil, sig_H_OH_sil, eps_H_OH_sil, q_H_OH_sil, &
                         sig_O_H2O, eps_O_H2O, q_O_H2O, sig_H_H2O, eps_H_H2O, q_H_H2O, &
                         sig_O_CO2, eps_O_CO2, q_O_CO2, sig_C_CO2, eps_C_CO2, q_C_CO2, &
                         sig_O_O2, eps_O_O2, q_O_O2, sig_N_N2, eps_N_N2, q_N_N2, &
                         bondl_SiO,bondl_OH,bondl_O2, bondl_N2, bondl_CO2,angle_SiOH
   USE dreiding_parameters_mod
   USE charmm_mod
   implicit none
   real(wp),parameter:: sigLJ_OH = sig_O_OH_sil
   real(wp),parameter:: angle_SiOSi = (140.0_wp/180.0_wp)*pi
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/Angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
! United atom H2O
   real(wp),parameter:: sigLJ_H2O = 3.166_wp/Angstrom
! United atom Si(OH)4
   real(wp),parameter:: r_tet = bondl_SiO + sigLJ_OH*0.5_wp
!
   real(wp),parameter:: q_Ox_sil = q_O_OH_sil
! LJ parameters
   integer,parameter:: ntyplj = ntyp  ! number of LJ types
   real(wp),parameter:: sigLJ(0:ntyplj) = (/ 0.0_wp, &
                        sig_Si_sil, &
                        sig_O_sil, &
                        sig_O_OH_sil, &
                        sig_H_OH_sil, &
                        sig_O_H2O, &
                        sig_H_H2O, &
                        sig_O_CO2, &
                        sig_C_CO2, &
                        sig_N_N2, &
                        sig_O_O2, &
                        sig_C2_sil, &
                        sig_C3_sil, &
                        sig_H_HC, &
                        sig_O_OC_Sil,&
                        sig_Ox_sil, &
                        sig_Ch2_sil, &
                        sig_Ch3_sil /)
   real(wp),parameter:: sigLJ_2(0:ntyplj) = sigLJ*0.5_wp
!
   real(wp),parameter:: epsLJ(0:ntyplj) = (/ 0.0_wp, &
                        eps_Si_sil, &
                        eps_O_sil, &
                        eps_O_OH_sil, &
                        eps_H_OH_sil, &
                        eps_O_H2O, &
                        eps_H_H2O, &
                        eps_O_CO2, &
                        eps_C_CO2, &
                        eps_N_N2, &
                        eps_O_O2, &
                        eps_C2_sil, &
                        eps_C3_sil, &
                        eps_H_HC, &
                        eps_O_OC_sil, &
                        eps_Ox_sil, &
                        eps_Ch2_sil, &
                        eps_Ch3_sil /)
!  integer,parameter:: in_N2charge = 16
!  integer,parameter:: iO_O2charge = 17
   real(wp),parameter:: qi(0:20) = (/ 0.0_wp, &
                        q_Si_sil, &
                        q_O_sil, &
                        q_O_OH_sil, &
                        q_H_OH_sil, &
                        q_O_H2O, &
                        q_H_H2O, &
                        q_O_CO2, &
                        q_C_CO2, &
                        q_N_N2, &
                        q_O_O2, &
                        q_C2_sil, &
                        q_C3_sil, &
                        q_H_HC, & 
                        q_O_OC_sil, &
                        q_Ox_sil, &
                        q_Ch2_sil, &
                        q_Ch3_sil, &
                        q_OH_sil, &
                -2.0_wp*q_N_N2, &
                -2.0_wp*q_O_O2 /)
END MODULE constants_mod

!!>include 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,Angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:)
   real(wp),allocatable:: Ox_xyz(:,:,:),r_O2_uc(:,:,:),dr_O2(:,:,:),rcm_O2_uc(:,:),vcm_O2(:,:)
   real(wp),allocatable:: Ox_fxyz(:,:,:),Ox_vxyz(:,:,:),image_Ox_xyz(:,:,:)
   real(wp),allocatable:: N2_xyz(:,:,:),r_N2_uc(:,:,:),dr_N2(:,:,:),rcm_N2_uc(:,:),vcm_N2(:,:)
   real(wp),allocatable:: N2_fxyz(:,:,:),N2_vxyz(:,:,:),image_N2_xyz(:,:,:)
   real(wp),allocatable:: CO2_xyz(:,:,:),r_CO2_uc(:,:,:),dr_CO2(:,:,:),rcm_CO2_uc(:,:),vcm_CO2(:,:)
   real(wp),allocatable:: CO2_fxyz(:,:,:),CO2_vxyz(:,:,:),image_CO2_xyz(:,:,:)
   real(wp),allocatable:: GAS_xyz(:,:,:),r_GAS_uc(:,:,:),dr_GAS(:,:,:),rcm_GAS_uc(:,:),vcm_GAS(:,:)
   real(wp),allocatable:: GAS_fxyz(:,:,:),GAS_vxyz(:,:,:),image_GAS_xyz(:,:,:)
   real(wp),allocatable:: h2o_xyz(:,:,:),r_h2o_uc(:,:,:),dr_h2o(:,:,:),rcm_h2o_uc(:,:),vcm_h2o(:,:)
   real(wp),allocatable:: h2o_fxyz(:,:,:),h2o_vxyz(:,:,:),image_h2o_xyz(:,:,:)
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
      write(iu,'(a6,1x,3(1x,f12.4))') 'REMARK',boxl*Angstrom
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


!!>include 'connectivity_mod.f90'

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
   !
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


   PURE FUNCTION OH_groups_OK(iOH1,iOH2)
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

    PURE FUNCTION is_teos_inter_group(iSi)
!
!     Is the Silicon iSi (attached to a bridging Oxygen)
!     in a 'teos_intermediate' group? i.e.
!                                     \
!     is it bonded to 3(OC)       OC--Si--OC
!     groups?                          |
!                                      OC
!
!                                     \
!     is it bonded to 2(OC),      OC--Si--OH
!     1(OH) groups?                     |
!                                      OC
!
!                                     \
!     is it bonded to 1(OC),      OC--Si--OH
!     2(OH) groups?                     |
!                                      OH
!
      integer,intent(in):: iSi
      logical:: is_teos_inter_group
      integer:: k,j,nOc1,nOH1
      is_teos_inter_group = .FALSE.
      nOc1 = 0
      nOH1 = 0
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         select case (atom(j))
         case(iOxygen)
            cycle
         case(iOxygenC)
            nOc1 = nOc1 + 1
         case(iOxygenH)
            nOH1 = nOH1 + 1
         end select
      end do
      if ((nOc1 == 3).or.(nOc1 == 2 .and. nOH1 == 1).or. &
          (nOc1 == 1 .and. nOH1 == 2)) then
         is_teos_inter_group = .TRUE.
      end if
   END FUNCTION is_teos_inter_group


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
            call write_atom_info(6,iSi)
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

  SUBROUTINE delete_group_C2C3(iSi,iOc)
!     Delete Ch2 and Ch3 united atoms 
      USE sort_mod,only:shell
      USE atom_types_mod
      integer,intent(inout):: iSi,iOc
      integer:: j,i,ii,jj
      integer:: ilst(7),ni
      ni = 0
      ilst = 0
      loop_1: do i = 1,ncmax(atom(iOc))
         ii = proximity(iOc,i)
         if (ii == iSi) cycle loop_1
         if (ii == 0) STOP 'error: delete_group_Ch2-Ch3'
         ni = ni + 1
         ilst(ni) = ii
         loop_2: do j = 1,ncmax(atom(ii))
            jj = proximity(ii,j)
            if (jj == iOc) cycle loop_2
            ni = ni + 1
            ilst(ni) = jj
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
   END SUBROUTINE delete_group_C2C3



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

   SUBROUTINE delete_OxygenH(iSi,iOH)
      USE sort_mod, only:shell
      implicit none
      integer,intent(inout):: iSi
      integer,intent(in)::iOH
      integer:: i,ii
      ! delete only H  attached with OH
      do i = 1,ncmax(atom(iSi))
         ii = proximity(iSi,i)
         if (ii == 0) cycle
         if (ii == iOH) exit
      end do
!      if (i > ncmax(atom(iOH))) STOP 'delete_Hydrogen: error'
      if (iSi == natom) then ! update the silicon number
         iSi = ii
      end if
      call delete_atom(ii)
   END SUBROUTINE delete_OxygenH

   SUBROUTINE del_teos_interm(iSi)
!Delete teos intermediate: Si(OC2H5)3(OH),Si(OC2H5)2(OH)2,Si(OC2H5)(OH)3 
      USE sort_mod,only:shell
      USE atom_types_mod
      integer,intent(in):: iSi
      integer:: k,j,i,l,ii,jj,kk,ll
      integer:: ilst(26),ni
      ni = 0
      ilst = 0
      loop_1: do i = 1,ncmax(atom(iSi))
         ii = proximity(iSi,i)         
         if (ii == 0) STOP 'error: del_tin_oc3oh'
         ni = ni + 1
         ilst(ni) = ii
         print*,'ii',ilst(ni)
         loop_2: do j = 1,ncmax(atom(ii))
            jj = proximity(ii,j)
            if (jj == 0)cycle loop_2
            if (jj == iSi) cycle loop_2
            ni = ni + 1
            ilst(ni) = jj
            print*,'jj',ilst(ni)
            loop_3: do k = 1,ncmax(atom(jj))
               kk = proximity(jj,k)
               if (kk == ii) cycle loop_3
               ni = ni + 1
               ilst(ni) = kk 
               print*,'kk',ilst(ni)           
               loop_4: do l = 1,ncmax(atom(kk))
                  ll = proximity(kk,l)
                  if (ll == jj) cycle loop_4
                  ni = ni + 1
                  ilst(ni) = ll
                  print*,'ll',ilst(ni)
               end do loop_4
            end do loop_3 
         end do loop_2
      end do loop_1
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE del_teos_interm

 SUBROUTINE del_teos_interm2(iSi)
!Delete teos intermediate: Si(OCh2Ch3)3(OH),Si(OCh2Ch3)2(OH)2,Si(OCh2Ch3)(OH)3 
      USE sort_mod,only:shell
      USE atom_types_mod
      integer,intent(in):: iSi
      integer:: k,j,i,ii,jj,kk
      integer:: ilst(11),ni
      ni = 0
      ilst = 0
      loop_1: do i = 1,ncmax(atom(iSi))
         ii = proximity(iSi,i)         
         if (ii == 0) STOP 'error: del_tin_oc3oh'
         ni = ni + 1
         ilst(ni) = ii
         print*,'ii',ilst(ni)
         loop_2: do j = 1,ncmax(atom(ii))
            jj = proximity(ii,j)
            if (jj == 0)cycle loop_2
            if (jj == iSi) cycle loop_2
            ni = ni + 1
            ilst(ni) = jj
            print*,'jj',ilst(ni)
            loop_3: do k = 1,ncmax(atom(jj))
               kk = proximity(jj,k)
               if (kk == ii) cycle loop_3
               ni = ni + 1
               ilst(ni) = kk 
               print*,'kk',ilst(ni)           
            end do loop_3 
         end do loop_2
      end do loop_1
      ni = ni + 1
      ilst(ni) = iSi
      call shell(ni,ilst)
      do i = ni,1,-1
         call delete_atom(ilst(i))
      end do
   END SUBROUTINE del_teos_interm2

       
   SUBROUTINE delete_groups(iSi)
      USE atom_types_mod
      integer,intent(in)::iSi
      integer:: i,j,nOc1,nOH1
      nOc1 = 0
      nOH1 = 0
      do i = 1,ncmax(atom(iSi))
         j = proximity(iSi,i)
         select case (atom(j)) 
         case(iOxygenC)
            nOc1 = nOc1 + 1
         case(iOxygenH)
            nOH1 = nOH1 + 1
         case default 
            stop'unkown molecule found'
         end select
      end do
      if ((nOc1 == 3.and.nOH1 ==1).or.(nOc1 == 2 .and. nOH1 == 2).or. &
          (nOc1 == 1 .and. nOH1 == 3)) then
          call del_teos_interm2(iSi)
      end if 
      if (nOH1 == 4)then
        call delete_group(iSi)
      end if       

   END SUBROUTINE delete_groups



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


   SUBROUTINE mt1_lst(iSi,itype,n,Lst)
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
    END SUBROUTINE mt1_lst


    SUBROUTINE mt2_lst(iSi,itype1,itype2,n1,Lst1,n2,Lst2)
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
   END SUBROUTINE mt2_lst


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

!!>include 'check_structure_mod_TEOS.f90'

MODULE check_structure_mod
   USE global_vars_mod
   USE atom_types_mod
   implicit none
CONTAINS

SUBROUTINE check_proximity2(ifirst,ilast,atoms,conect,consistent,ierr)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: atoms(:),conect(:,:)
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k,m,n,ia,a1,a2,i1,i2
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         ia = atoms(i)
         select case(ia)
         case(iHydrogenC)  ! Check H is bonded to 1 Carbon2 or Carbon3
            m = conect(i,1)
            if (m == 0) GOTO 100
            if (atoms(m) == iSilicon) GOTO 100
            if (atoms(m) == iOxygenC) GOTO 100
            if (atoms(m) == iOxygenF) GOTO 100
            if (atoms(m) == iOxygen) GOTO 100
            if (atoms(m) == iOxygenH) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
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
         case(iOxygenC)  ! OxygenC must be bonded to Si and Carbon2
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            i2 = conect(i,2)
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (.NOT.((a1 == iSilicon .and. a2 == iCarbonH2).or. &
                      (a2 == iSilicon .and. a1 == iCarbonH2))) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
         case(iOxygenF) !Check free oxygen
            n = conect(i,1)
            if (atoms(n) /= iSilicon) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
         case(iCarbon2)  ! Check C2 is bonded to 2 HC, OC and C3 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iOxygen) GOTO 100
               if (atoms(k) == iOxygenH) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100
               if (atoms(k) == iOxygenF) GOTO 100
            end do
         case(iCarbon3)  ! Check C3 is bonded to 3 HC and C2 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iOxygen) GOTO 100
               if (atoms(k) == iOxygenH) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100
               if (atoms(k) == iOxygenF) GOTO 100
               if (atoms(k) == iOxygenC) GOTO 100
            end do

         case(iCarbonH2)  ! Check Ch2 is bonded to Oc and Ch3 
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            i2 = conect(i,2)
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (.NOT.((a1 == iOxygenC .and. a2 == iCarbonH3).or. &
                      (a2 == iOxygenC .and. a1 == iCarbonH3))) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
            
         case(iCarbonH3)  ! Check C3 is bonded to Ch2 atoms
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            if (a1 /= iCarbonH2)GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100    
            
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
            STOP 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ierr = i
   END SUBROUTINE check_proximity2


   SUBROUTINE check_proximity2_H(ifirst,ilast,atoms,conect,consistent,ierr)
!     Check that atoms are bonded to the correct types of atoms.
!     e.g. NO Si-Si bonds, etc.
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: atoms(:),conect(:,:)
      logical,intent(out):: consistent
      integer,intent(out):: ierr
      integer:: i,j,k,m,n,ia,a1,a2,i1,i2
      consistent = .true.
      ierr = 0
      do i = ifirst,ilast
         ia = atoms(i)
         select case(ia)
         case(iHydrogen)! Check H is bonded to 1 OH Oxygen
            k = conect(i,1)
            if (k == 0) GOTO 100
            if (atoms(k) /= iOxygenH) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
         case(iHydrogenC)  ! Check H is bonded to 1 Carbon2 or Carbon3
            m = conect(i,1)
            if (m == 0) GOTO 100
            if (atoms(m) == iHydrogen) GOTO 100
            if (atoms(m) == iSilicon) GOTO 100
            if (atoms(m) == iOxygenC) GOTO 100
            if (atoms(m) == iOxygenF) GOTO 100
            if (atoms(m) == iOxygen) GOTO 100
            if (atoms(m) == iOxygenH) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
         case(iSilicon)  ! Check Si is bonded to 4 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100   ! but not Si
               if (atoms(k) == iHydrogen) GOTO 100
            end do
         case(iOxygen)  ! Bridging Oxygen must be bonded to Si
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            if (a1 /= iSilicon) GOTO 100
!if (i > nseed) then  ! non seed bridging O
            i2 = conect(i,2)  !  must be bonded to  2 Silicons
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (a2 /= iSilicon) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
!end if
         case(iOxygenC)  ! OxygenC must be bonded to Si and Carbon2
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
            i2 = conect(i,2)
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (.NOT.((a1 == iSilicon .and. a2 == iCarbon2).or. &
                      (a2 == iSilicon .and. a1 == iCarbon2))) GOTO 100
            if (any(conect(i,3:) /= 0)) GOTO 100
         case(iOxygenF) !Check free oxygen
            n = conect(i,1)
            if (atoms(n) /= iSilicon) GOTO 100
            if (any(conect(i,2:) /= 0)) GOTO 100
         case(iCarbon2)  ! Check C2 is bonded to 2 HC, OC and C3 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iHydrogen) GOTO 100
               if (atoms(k) == iOxygen) GOTO 100
               if (atoms(k) == iOxygenH) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100
               if (atoms(k) == iOxygenF) GOTO 100
            end do
         case(iCarbon3)  ! Check C3 is bonded to 3 HC and C2 atoms
            do j = 1,ncmax(ia)
               k = conect(i,j)
               if (k == 0) GOTO 100
               if (atoms(k) == iHydrogen) GOTO 100
               if (atoms(k) == iOxygen) GOTO 100
               if (atoms(k) == iOxygenH) GOTO 100
               if (atoms(k) == iSilicon) GOTO 100
               if (atoms(k) == iOxygenF) GOTO 100
               if (atoms(k) == iOxygenC) GOTO 100
            end do
         case(iOxygenH)  ! OH Oxygen must be bonded to Si and H
            i1 = conect(i,1)
            if (i1 == 0) GOTO 100
            a1 = atoms(conect(i,1))
!if (i > nseed) then
            i2 = conect(i,2)
            if (i2 == 0) GOTO 100
            a2 = atoms(conect(i,2))
            if (i2 == 0) GOTO 100
            if (.NOT.((a1 == iHydrogen .and. a2 == iSilicon).or. &
                      (a2 == iHydrogen .and. a1 == iSilicon))) GOTO 100
            if ( any(conect(i,3:) /= 0) ) GOTO 100
!end if
         case default
            STOP 'unknown atom type'
         end select
      end do
      RETURN
100   consistent = .false.
      ierr = i
   END SUBROUTINE check_proximity2_H

END MODULE check_structure_mod

!!>include 'bond_angle_types_mod_TEOS.f90'

MODULE bond_angle_types_mod
    USE precision_mod
    USE atom_types_mod
    USE constants_mod
    implicit none

CONTAINS

   SUBROUTINE bond_type(ia1,ia2,abond,kbond)
      integer,intent(in):: ia1,ia2
      real(wp),intent(out):: kbond,abond
      if (compare_bond(ia1,ia2,iSilicon,iOxygen)) then
         kbond = kSiO
         abond = aSiO
      else if (compare_bond(ia1,ia2,iSilicon,iOxygenH)) then
         kbond = kSiO
         abond = aSiO
      else if(compare_bond(ia1,ia2,iSilicon,iSilicon))then
         kbond = KSiSi
         abond = ASiSi
         STOP 'Si-Si bond not allowed'
      else if(compare_bond(ia1,ia2,iSilicon,iOxygenC))then
         kbond = KSiOc
         abond = ASiOc
      else if(compare_bond(ia1,ia2,iCarbon2,iCarbon3))then
         kbond = KC2C3
         abond = AC2C3
      else if(compare_bond(ia1,ia2,iOxygenC,iCarbon2))then
         kbond = KOcC2
         abond = AOcC2
      else if(compare_bond(ia1,ia2,iCarbon2,iHydrogenC))then
         kbond = KC2Hc
         abond = AC2Hc
      else if(compare_bond(ia1,ia2,iCarbon3,iHydrogenC))then
         kbond = KC3Hc
         abond = AC3Hc
      else if(compare_bond(ia1,ia2,iCarbonH2,iCarbonH3))then
         kbond = KC2C3
         abond = AC2C3
      else if(compare_bond(ia1,ia2,iOxygenC,iCarbonH2))then
         kbond = KOcC2
         abond = AOcC2
            else
         STOP 'unknown bond'
      end if
   END SUBROUTINE bond_type


   SUBROUTINE angle_type(ia1,ia2,ia3,ctheta,kang)
      integer,intent(in):: ia1,ia2,ia3
      real(wp),intent(out):: ctheta,kang
      if ( compare_angle(ia1,ia2,ia3,iOxygen,iSilicon,iOxygen) .or.&
           compare_angle(ia1,ia2,ia3,iOxygen,iSilicon,iOxygenH) .or.&
           compare_angle(ia1,ia2,ia3,iOxygen,iSilicon,iOxygenC) .or.&
           compare_angle(ia1,ia2,ia3,iOxygenH,iSilicon,iOxygenH) .or.&
           compare_angle(ia1,ia2,ia3,iOxygenC,iSilicon,iOxygenH) .or.&
           compare_angle(ia1,ia2,ia3,iOxygenC,iSilicon,iOxygenC) ) then
         kang = KOSiO
         ctheta = cos_OSiO
      else if (compare_angle(ia1,ia2,ia3,iSilicon,iOxygen,iSilicon)) then
         kang = KSiOSi
         ctheta = cos_SiOSi
      else if(compare_angle(ia1,ia2,ia3,iSilicon,iOxygenC,iCarbon2)) then
         kang = kSiOcC2
         ctheta = cos_SiOC
      else if(compare_angle(ia1,ia2,ia3,iOxygenC,iCarbon2,iCarbon3)) then
         kang = KOC2C3
         ctheta = cos_OC2C3
      else if(compare_angle(ia1,ia2,ia3,iOxygenC,iCarbon2,iHydrogenC)) then
         kang = KOcC2Hc
         ctheta = cos_OSiO
      else if(compare_angle(ia1,ia2,ia3,iCarbon2,iCarbon3,iHydrogenC)) then
         kang = KC2C3Hc
         ctheta = cos_OSiO
      else if(compare_angle(ia1,ia2,ia3,iCarbon3,iCarbon2,iHydrogenC)) then
         kang = KHcC3Hc
         ctheta = cos_OSiO
      else if(compare_angle(ia1,ia2,ia3,iOxygenC,iCarbonH2,iCarbonH3)) then
         kang = KOC2C3
         ctheta = cos_OC2C3
      else if(compare_angle(ia1,ia2,ia3,iSilicon,iOxygenC,iCarbonH2)) then
         kang = kSiOcC2
         ctheta = cos_SiOC
      else if(compare_angle(ia1,ia2,ia3,iHydrogenC,iCarbon2,iHydrogenC).or.&
              compare_angle(ia1,ia2,ia3,iHydrogenC,iCarbon3,iHydrogenC) ) then
         kang = KHCH
         ctheta = cos_HCH
      else
         print *,ia1,ia2,ia3
         STOP 'unknown angle'
      end if
   END SUBROUTINE angle_type


   logical pure function compare_bond(i1,i2,typ1,typ2)
      integer,intent(in):: i1,i2,typ1,typ2
      compare_bond = (i1==typ1.and.i2==typ2).or. &
                     (i1==typ2.and.i2==typ1)
   end function compare_bond


   logical pure function compare_angle(i1,i2,i3,typ1,typ2,typ3)
      integer,intent(in):: i1,i2,i3,typ1,typ2,typ3
      compare_angle = (i1==typ1.and.i2==typ2.and.i3==typ3).or. &
                      (i1==typ3.and.i2==typ2.and.i3==typ1)
   end function compare_angle

END MODULE bond_angle_types_mod


!!>include 'bond_angle_list_mod.f90'

MODULE bond_angle_list_mod
   USE precision_mod
   USE global_vars_mod
   USE bond_angle_types_mod
   USE connectivity_mod, only: proximity
   USE atom_types_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   real(wp),allocatable:: kbond(:),abond(:)
   integer:: nbondtot,nang
!
   TYPE bond_angle_list
      integer:: n
      integer:: nbond, nang
      integer,allocatable:: ibond(:,:),iang(:,:)
      real(wp),allocatable:: kbond(:),abond(:)
      real(wp),allocatable:: kang(:),cang(:)
   END TYPE bond_angle_list
!
CONTAINS

   SUBROUTINE new_bond_angle_list(T,n)
      type(bond_angle_list),intent(out):: T
      integer,intent(in):: n
      T%n = n
      allocate( T%ibond(2,n*4) )
      allocate( T%iang(3,n*6) )
      allocate( T%kbond(n*4),T%abond(n*4) )
      allocate( T%kang(n*6),T%cang(n*6) )
      T%nbond = 0
      T%nang = 0
   END SUBROUTINE new_bond_angle_list


   SUBROUTINE bond_list()
      integer:: ia1,ia2,j
      nbondtot = 0
      if (.not.allocated(ibond)) allocate( ibond(2,natom_max*4))
      if (.not.allocated(nbond)) allocate( nbond(natom_max) )
!     if (.not.allocated(kbond)) allocate( kbond(natom_max),abond(natom_max) )
!
      do ia1 = 1,natom
         nbond(ia1) = 0
         do j = 1,ncmax(atom(ia1))
            ia2 = proximity(ia1,j)
            if (ia2 > 0) nbond(ia1) = nbond(ia1) + 1
            if (ia2 > ia1) then
               nbondtot = nbondtot + 1
               ibond(1,nbondtot) = ia1
               ibond(2,nbondtot) = ia2
               !call bond_type(atom(ia1),atom(ia2),abond(nbondtot),kbond(nbondtot))
            end if
         end do
      end do
   END SUBROUTINE bond_list


   SUBROUTINE set_bond_angle_lists()
      integer:: ia1,ia2,ia3,ib,ic
!
      if (.not.allocated(ibond)) then
         allocate( ibond(2,natom_max*4) )
         allocate( iang(3,natom_max*6) )
         allocate( kbond(natom_max*4),abond(natom_max*4))
         allocate( kang(natom_max*6),ctheta(natom_max*6))
      end if
      nbondtot = 0
      nang = 0
      do ia1 = 1,natom
      do ib = 1,ncmax(atom(ia1))
         ia2 = proximity(ia1,ib)
         if (ia2 == 0) cycle
         if (ia1 < ia2) then
            nbondtot = nbondtot + 1
            ibond(1,nbondtot) = ia1
            ibond(2,nbondtot) = ia2
            call bond_type(atom(ia1),atom(ia2),abond(nbondtot),kbond(nbondtot))
            ! print *,ia1,ia2
            ! print *,kbond(nbondtot),abond(nbondtot)
         end if
         do ic = 1,ncmax(atom(ia2))
            ia3 = proximity(ia2,ic)
            if (ia3 == 0) cycle
            if (ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
            if (ia1 < ia3) then
               nang = nang + 1
               iang(1,nang) = ia1
               iang(2,nang) = ia2
               iang(3,nang) = ia3
               call angle_type(atom(ia1),atom(ia2),atom(ia3),ctheta(nang),kang(nang))
            end if
         end do
      end do
      end do
!      do i = 1, nang
!         if (atom(iang(2,i)) == iSilicon) then
!            kang(i) = KOSiO
!            ctheta(i) = cos_OSiO
!         else
!            kang(i) = KSiOSi
!            ctheta(i) = cos_SiOSi
!         end if
!      end do
   END SUBROUTINE set_bond_angle_lists


   SUBROUTINE set_bond_angle_lists2(ifirst,ilast,conect,atype,nbond,ibond,kbond,abond,nang,iang,kang,cang)
      integer,intent(in):: conect(:,:),atype(:),ifirst,ilast
      integer,intent(out):: nbond,ibond(:,:),nang,iang(:,:)
      real(wp),intent(out):: kbond(:),abond(:),kang(:),cang(:)
      integer:: ia1,ia2,ia3,ib,ic
!
      nbond = 0
      nang = 0
      do ia1 = ifirst,ilast
      do ib = 1,ncmax(atype(ia1))
         ia2 = conect(ia1,ib)
         if (ia2 == 0) cycle
         if (ia1 < ia2) then
            nbond = nbond + 1
            ibond(1,nbond) = ia1
            ibond(2,nbond) = ia2
            call bond_type(atype(ia1),atype(ia2),abond(nbond),kbond(nbond))
            ! print *,ia1,ia2
            ! print *,kbond(nbond),abond(nbond)
         end if
         do ic = 1,ncmax(atype(ia2))
            ia3 = conect(ia2,ic)
            if (ia3 == 0) cycle
            if (ia1 == ia3) cycle   ! if atom_1 == atom_3 skip ang
            if (ia1 < ia3) then
               nang = nang + 1
               iang(1,nang) = ia1
               iang(2,nang) = ia2
               iang(3,nang) = ia3
               call angle_type(atype(ia1),atype(ia2),atype(ia3),cang(nang),kang(nang))
            end if
         end do
      end do
      end do
   END SUBROUTINE set_bond_angle_lists2


   SUBROUTINE atom_bond_angle(ia,nbnd,ibnd,kbnd,abnd,nan,ian,kan,can)
      integer,intent(in):: ia
      integer,intent(out):: nbnd,ibnd(:,:),nan,ian(:,:)
      real(wp),intent(out):: kbnd(:),abnd(:),kan(:),can(:)
      integer:: jj,kk,j,k
!
! the angles with ia as the centre atom: j--ia--k
      nan = 0
      do jj = 1, ncmax(atom(ia))-1
         j = proximity(ia,jj)
         if (j == 0) cycle
         do kk = jj+1, ncmax(atom(ia))
            k = proximity(ia,kk)
            if (k == 0) cycle
            nan = nan + 1
            ian(1,nan) = j
            ian(2,nan) = ia
            ian(3,nan) = k
            call angle_type(atom(j),atom(ia),atom(k),can(nan),kan(nan))
         end do
      end do
! the bonds involving atom ia and
! the angles with ia as an end atom: ia--j--k
      nbnd = 0
      do jj = 1, ncmax(atom(ia))
         j = proximity(ia,jj)
         if (j == 0) cycle
         nbnd = nbnd + 1
         ibnd(1,nbnd) = ia
         ibnd(2,nbnd) = j
         call bond_type(atom(ia),atom(j),abnd(nbnd),kbnd(nbnd))
         do kk = 1, ncmax(atom(j))
            k = proximity(j,kk)
            if (k==ia .or. k==0) cycle
            nan = nan + 1
            ian(1,nan) = ia
            ian(2,nan) = j
            ian(3,nan) = k
            call angle_type(atom(ia),atom(j),atom(k),can(nan),kan(nan))
         end do
      end do
   END SUBROUTINE atom_bond_angle


   SUBROUTINE list_bond_angle(nl,list,nbnd,ibnd,kbnd,abnd,nan,ian,kan,can)
!
! Given a list of atoms create the bond & angle lists for the interactions
! of the list atoms with each other & with the 'boundary' atoms to the list.
!
      integer,intent(in):: nl,list(:)
      integer,intent(out):: nbnd,ibnd(:,:),nan,ian(:,:)
      real(wp),intent(out):: kbnd(:),abnd(:),kan(:),can(:)
      integer:: i,jj,kk,j,k,il,indj,indk
!
      nbnd = 0
      nan = 0
      list_loop: do il = 1, nl
      i = list(il)
!
! the angles with i as the centre atom: j--i--k
      do jj = 1, ncmax(atom(i))-1
         j = proximity(i,jj)
         if (j == 0) cycle
         do kk = jj+1, ncmax(atom(i))
            k = proximity(i,kk)
            if (k == 0) cycle
            nan = nan + 1
            ian(1,nan) = j
            ian(2,nan) = i
            ian(3,nan) = k
            call angle_type(atom(j),atom(i),atom(k),can(nan),kan(nan))
         end do
      end do
! the bonds involving atom i and
! the angles with i as an end atom: i--j--k where j is not in list
      do jj = 1, ncmax(atom(i))
         j = proximity(i,jj)
         if (j == 0) cycle
         indj = indx(nl,list,j)
         if (indj > 0 .and. indj < il) cycle
         nbnd = nbnd + 1
         ibnd(1,nbnd) = i
         ibnd(2,nbnd) = j
         call bond_type(atom(i),atom(j),abnd(nbnd),kbnd(nbnd))

         ! if (ANY(list(1:nl) == j)) cycle
         ! if (in_list(j,list,nl)) cycle
         if (indj > 0) cycle
         !
         ! if j is not in list
         !          |     k
         !  list    |   /
         !       i--|--j----k
         !          |   \
         !          |    k  outside list
         !
         do kk = 1, ncmax(atom(j))
            k = proximity(j,kk)
            if (k==i .or. k==0) cycle
            indk = indx(nl,list,k)
            if (indk > 0) then
               !
               !  list   \
               !       i--\--j
               !           \/   j  outside list
               !           /\   BUT k inside list
               !          k  \
               if (indk < il) cycle
            end if
            nan = nan + 1
            ian(1,nan) = i
            ian(2,nan) = j
            ian(3,nan) = k
            call angle_type(atom(i),atom(j),atom(k),can(nan),kan(nan))
         end do
      end do
!
      end do list_loop
   contains
      pure function indx(nl,list,k)
         integer:: indx
         integer,intent(in):: nl,list(:),k
         integer:: i
         indx = 0
         do i = 1,nl
            if (list(i) == k) then
               indx = i
               RETURN
            end if
         end do
      end function indx
   END SUBROUTINE list_bond_angle


   SUBROUTINE list_bond_angle_sort(nl,lst,nbnd,ibnd,kbnd,abnd,nan,ian,kan,can)
      USE sort_mod, only: qsort
!
! Given a list of atoms create the bond & angle lists for the interactions
! of the list atoms with each other & with the 'boundary' atoms to the list.
!
      integer,intent(in):: nl,lst(:)
      integer,intent(out):: nbnd,ibnd(:,:),nan,ian(:,:)
      real(wp),intent(out):: kbnd(:),abnd(:),kan(:),can(:)
      integer:: list(nl)
      integer:: i,jj,kk,j,k,il,indj,indk
!
      call qsort(nl,lst(1:nl),list)
      nbnd = 0
      nan = 0
      list_loop: do il = 1, nl
      i = list(il)
!
! the angles with i as the centre atom: j--i--k
      do jj = 1, ncmax(atom(i))-1
         j = proximity(i,jj)
         if (j == 0) cycle
         do kk = jj+1, ncmax(atom(i))
            k = proximity(i,kk)
            if (k == 0) cycle
            nan = nan + 1
            ian(1,nan) = j
            ian(2,nan) = i
            ian(3,nan) = k
            call angle_type(atom(j),atom(i),atom(k),can(nan),kan(nan))
         end do
      end do
! the bonds involving atom i and
! the angles with i as an end atom: i--j--k where j is not in list
      do jj = 1, ncmax(atom(i))
         j = proximity(i,jj)
         if (j == 0) cycle
         indj = indx(nl,list,j)
         if (indj > 0 .and. indj < il) cycle
         nbnd = nbnd + 1
         ibnd(1,nbnd) = i
         ibnd(2,nbnd) = j
         call bond_type(atom(i),atom(j),abnd(nbnd),kbnd(nbnd))

         ! if (ANY(list(1:nl) == j)) cycle
         ! if (in_list(j,list,nl)) cycle
         if (indj > 0) cycle
         !
         ! if j is not in list
         !          |     k
         !  list    |   /
         !       i--|--j----k
         !          |   \
         !          |    k  outside list
         !
         do kk = 1, ncmax(atom(j))
            k = proximity(j,kk)
            if (k==i .or. k==0) cycle
            indk = indx(nl,list,k)
            if (indk > 0) then
               !
               !  list   \
               !       i--\--j
               !           \/   j  outside list
               !           /\   BUT k inside list
               !          k  \
               if (indk < il) cycle
            end if
            nan = nan + 1
            ian(1,nan) = i
            ian(2,nan) = j
            ian(3,nan) = k
            call angle_type(atom(i),atom(j),atom(k),can(nan),kan(nan))
         end do
      end do
!
      end do list_loop
   contains
      pure function indx(nl,list,k)
         integer:: indx
         integer,intent(in):: nl,list(:),k
         integer:: i
         indx = 0
         if (k < list(1) .or. k > list(nl)) RETURN
         do i = 1,nl
            if (list(i) == k) then
               indx = i
               RETURN
            end if
         end do
      end function indx
   END SUBROUTINE list_bond_angle_sort

END MODULE bond_angle_list_mod

!!>include 'box_frames_mod.f90'

MODULE BOX_FRAMES_MOD
   USE precision_mod, only: wp
   USE coordinates_mod
   implicit none
   integer,public:: imve, iframe
   public:: new_frame_file, write_box_frame
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
      write(unit=ctmp,fmt='(i5.5)') iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.xyz'
      inquire(file=trim(frame_file),iostat=ios,exist = lexists,opened = lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit=iu,file=trim(adjustl(frame_file)),form = 'formatted')
   END SUBROUTINE new_frame_file


   SUBROUTINE write_box_frame(iu,nmol,print_all)
      USE global_vars_mod
      USE atom_types_mod
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j,k
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 2*n_O2 +2*n_N2 +3*n_CO2 + n_atom_GAS*n_GAS + 2 + n_atom_h2o*n_h2o
            write(iu,'(a,6f14.8)')'Gases+Silica ',boxl*Angstrom
         else
            write(iu,*) 2*n_O2 +2*n_N2 +3*n_CO2 + n_atom_GAS*n_GAS + 2 + n_atom_h2o*n_h2o
            write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
         end if
      else
         write(iu,*) 2*n_O2 +2*n_N2 +3*n_CO2 + n_atom_GAS*n_GAS + 2 + n_atom_h2o*n_h2o
         write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
      end if
      do i = 1,n_O2
         write(iu,110)'O2 ',Ox_xyz(i,1,:)*Angstrom
         write(iu,110)'O2 ',Ox_xyz(i,2,:)*Angstrom
      end do
      do i = 1,n_N2
         write(iu,110)'N2 ',N2_xyz(i,1,:)*Angstrom
         write(iu,110)'N2 ',N2_xyz(i,2,:)*Angstrom
      end do
      do i = 1,n_CO2
         write(iu,110)'CO ',CO2_xyz(i,1,:)*Angstrom
         write(iu,110)'OC ',CO2_xyz(i,2,:)*Angstrom
         write(iu,110)'OC ',CO2_xyz(i,3,:)*Angstrom
      end do
      do i = 1,n_GAS
         do j = 1,n_atom_GAS
            write(iu,110) atom_name(gas_atom(i,j)),(GAS_xyz(i,j,k)*Angstrom,k = 1,NDIM)
         end do
      end do
      do i = 1,n_h2o
         do j = 1,n_atom_h2o
            write(iu,110) atom_name(h2o_atom(i,j)),(h2o_xyz(i,j,k)*Angstrom,k = 1,NDIM)
         end do
      end do

!
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,110) atom_name(atom(i)),(rxyz(i,j)*Angstrom,j = 1,NDIM)
         end do
         end if
      end if
      write(iu,110)'Ar ',(/-17.832,-17.832, -10.000/)
      write(iu,110)'Ar ',(/-17.832, -17.832, 100.000/)
!      write(iu,110)'Ar ',(/17.832, -17.832, 100.000/)
!      write(iu,110)'Ar ',  boxl2*Angstrom
!      write(iu,110)'Ar ', -boxl2*Angstrom
      write(iu,*)
      close(iu,status='keep')
110   format(a,6(1x,f12.6))
   END SUBROUTINE write_box_frame

   SUBROUTINE write_vel(iu,nmol,print_all)
      USE global_vars_mod
      USE atom_types_mod
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j,k
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 2*n_O2 +2*n_N2 +3*n_CO2 + n_atom_GAS*n_GAS + n_atom_h2o*n_h2o
            write(iu,'(a,6f14.8)')'Gases+Silica ',boxl*Angstrom
         else
            write(iu,*) 2*n_O2 +2*n_N2 +3*n_CO2 + n_atom_GAS*n_GAS + n_atom_h2o*n_h2o
            write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
         end if
      else
         write(iu,*) 2*nmol +2*nmol +3*nmol + n_atom_GAS*n_GAS + n_atom_h2o*n_h2o
         write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
      end if
      do i = 1,n_O2
         write(iu,*) (Ox_vxyz(i,1,k),k = 1,NDIM)
         write(iu,*) (Ox_vxyz(i,2,k),k = 1,NDIM)
      end do
      do i = 1,n_N2
         write(iu,*) (N2_vxyz(i,1,k),k = 1,NDIM)
         write(iu,*) (N2_vxyz(i,2,k),k = 1,NDIM)
      end do
      do i = 1,n_CO2
         write(iu,*) (CO2_vxyz(i,1,k),k = 1,NDIM)
         write(iu,*) (CO2_vxyz(i,2,k),k = 1,NDIM)
         write(iu,*) (CO2_vxyz(i,3,k),k = 1,NDIM)
      end do
      do i = 1,n_GAS
         do j = 1,n_atom_GAS
            write(iu,*) (GAS_vxyz(i,j,k),k = 1,NDIM)
         end do
      end do
      do i = 1,n_h2o
         do j = 1,n_atom_h2o
            write(iu,*) (h2o_vxyz(i,j,k),k = 1,NDIM)
         end do
      end do
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,*) (vxyz(i,k),k = 1,NDIM)
         end do
         end if
      end if
      write(iu,*)
      close(iu,status='keep')
   END SUBROUTINE write_vel

   SUBROUTINE read_vel(iu,nmol,print_all)
      USE global_vars_mod
      USE atom_types_mod
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      integer:: i,j,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp
      do i = 1,n_O2
         read(iu,*) Ox_vxyz(i,1,:)
         read(iu,*) Ox_vxyz(i,2,:)
      end do
      do i = 1,n_N2
         read(iu,*) N2_vxyz(i,1,:)
         read(iu,*) N2_vxyz(i,2,:)
      end do
      do i = 1,n_CO2
         read(iu,*) CO2_vxyz(i,1,:)
         read(iu,*) CO2_vxyz(i,2,:)
         read(iu,*) CO2_vxyz(i,3,:)
      end do
      do i = 1,n_GAS
         do j = 1,n_atom_GAS
            read(iu,*) GAS_vxyz(i,j,:)
         end do
      end do
      do i = 1,n_h2o
         do j = 1,n_atom_h2o
            read(iu,*) h2o_vxyz(i,j,:)
         end do
      end do
      if (present(print_all)) then
      if (print_all) then
         do i = 1,natom
            read(iu,*) vxyz(i,1:NDIM)
         end do
      end if
      end if
      !close(iu,status='keep')
   END SUBROUTINE read_vel


   SUBROUTINE read_box_frame(iu,nmol,print_all)
      USE global_vars_mod
      USE atom_types_mod
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      real(wp):: r(NDIM)
      integer:: i,j,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp
      do i = 1,n_O2
         read(iu,*) ctmp,r
         Ox_xyz(i,1,:) = r/Angstrom
         read(iu,*) ctmp,r
         Ox_xyz(i,2,:) = r/Angstrom
      end do
      do i = 1,n_N2
         read(iu,*) ctmp,r
         N2_xyz(i,1,:) = r/Angstrom
         read(iu,*) ctmp,r
         N2_xyz(i,2,:) = r/Angstrom
      end do
      do i = 1,n_CO2
         read(iu,*) ctmp,r
         CO2_xyz(i,1,:) = r/Angstrom
         read(iu,*) ctmp,r
         CO2_xyz(i,2,:) = r/Angstrom
         read(iu,*) ctmp,r
         CO2_xyz(i,3,:) = r/Angstrom
      end do
      do i = 1,n_GAS
         do j = 1,n_atom_GAS
            read(iu,*) ctmp,r
            GAS_xyz(i,j,:) = r/Angstrom
         end do
      end do
      do i = 1,n_h2o
         do j = 1,n_atom_h2o
            read(iu,*) ctmp,r
            h2o_xyz(i,j,:) = r/Angstrom
         end do
      end do
      if (present(print_all)) then
      if (print_all) then
         do i = 1,natom
            read(iu,*) ctmp,r
            rxyz(i,1:NDIM) = r/Angstrom
         end do
      end if
      end if
      !close(iu,status='keep')
   END SUBROUTINE read_box_frame

   SUBROUTINE write_frame(iu,ifirst,ilast)
      USE files_mod
      USE bond_angle_list_mod
      USE global_vars_mod, only: Angstrom
      USE connectivity_mod
      USE atom_types_mod, only: atom
      USE coordinates_mod
      integer,intent(in):: iu,ifirst,ilast
      integer:: i,ib,ia,iu2,nc
      logical:: lopen,fix_bonds
      character(len=132):: frame_file
      call WRITE_PDB_HEADER(iu,ilast-ifirst+1,boxl)
      call WRITE_PDB(iu,ifirst,ilast,atom,rxyz)
      call bond_list
!
      inquire(unit=iu,name=frame_file)
      nc = len_trim(frame_file)
      frame_file = frame_file(1:nc-3)//'scr'
      fix_bonds = .false.
      do i = 1,nbondtot
         ia = ibond(1,i)
         ib = ibond(2,i)
         if (ANY( ABS(rxyz(ia,1:3) - rxyz(ib,1:3)) > 10.0_wp/Angstrom )) then
            fix_bonds = .true.
            inquire(file=frame_file,opened=lopen)
            if (.not.lopen) then
               call get_free_file_unit(iu2)
               open(unit=iu2,file=trim(frame_file))
            end if
            write(iu2,'(a,i0,a,i0)')'select atomno=',ia,', atomno=',ib
            write(iu2,*)'wireframe off'
           ! write(iu2,'(a,i0,a,i0)')'delbond top ',ia-1,' ',ib-1
         end if
      end do
      if (fix_bonds) then
         write(iu2,*)
         close(iu2,status='keep')
      end if
      call write_conect_pdb(iu,ifirst,ilast,proximity)
      call write_box_pdb(iu,(/ -boxl2(1),-boxl2(2), 0.0_wp /)*Angstrom, &
                            (/  boxl2(1), boxl2(2), 0.5_wp /)*Angstrom,99990,'CR')
      write(iu,'(a/)')'END'
      close(iu,status='keep')
   END SUBROUTINE write_frame

END MODULE BOX_FRAMES_MOD

!!>include 'vverlet_mod.f90'

MODULE vverlet_mod
    USE precision_mod
    implicit none
CONTAINS

! First part of water/salt velocity Verlet algorithm
   SUBROUTINE vverlet_a(dt,massi,r,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: r(:,:),v(:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      r = r + dt*v + dtsq2*f*massi
      v = v + dt2*f*massi
   END SUBROUTINE vverlet_a

! Second part of velocity Verlet
   SUBROUTINE vverlet_b(dt,massi,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: v(:,:)
      v = v + dt*0.5_wp*f*massi
   END SUBROUTINE vverlet_b

   SUBROUTINE gas_vverlet_a(dt,massi,r,v,f,dr)
      real(wp),intent(in):: dt,massi(:,:,:),f(:,:,:)
      real(wp),intent(inout):: r(:,:,:),v(:,:,:)
      real(wp),intent(inout):: dr(:,:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      dr = dt*v + dtsq2*f*massi
      r = r + dr
      v = v + dt2*f*massi
   END SUBROUTINE gas_vverlet_a

   SUBROUTINE gas_vverlet_b(dt,massi,v,f)
      real(wp),intent(in):: dt,massi(:,:,:),f(:,:,:)
      real(wp),intent(inout):: v(:,:,:)
      v = v + dt*0.5_wp*f*massi
   END SUBROUTINE gas_vverlet_b

END MODULE vverlet_mod

!!>include 'nlist_mod.f90'

MODULE NLIST_MOD
   USE precision_mod,only:wp
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
 print*,'bottom_corner = ',bottom_corner 
      side = rupper - rlower
 print*,'side = ',side       
      ncelx = int(side(1)/Rc)
      ncely = int(side(2)/Rc)
      ncelz = int(side(3)/Rc)
      ncelt = ncelx*ncely*ncelz
print*,'ncelx =',ncelx 
print*,'ncely =',ncely 
print*,'ncelz =',ncelz 
print*,'ncelt =',ncelt 
      del_nl_x = side(1)/ncelx
      del_nl_y = side(2)/ncely
      del_nl_z = side(3)/ncelz      
      delmin = min(del_nl_x,del_nl_y,del_nl_z)
      delxi = 1.0_wp/del_nl_x
      delyi = 1.0_wp/del_nl_y
      delzi = 1.0_wp/del_nl_z
print*,'del_nl_x =',del_nl_x 
print*,'del_nl_y =',del_nl_y 
print*,'del_nl_z =',del_nl_z 
print*,'delmin =',delmin 
print*,'delxi =',delxi 
print*,'delyi =',delyi 
print*,'delzi =',delzi 
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
         !if (ic == 100)stop
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



!!>include 'verlet_list_nl.f90'

MODULE verlet_list_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom, natom_max
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

   SUBROUTINE INIT_VERLET_LIST(rv0,rcut0,nat_max)
      real(wp),intent(in):: rv0,rcut0
      integer,intent(in):: nat_max
      rv = rv0
      rcut = rcut0
      RV2 = RV*RV
      SKIN = RV - rcut
      SKIN2 = SKIN*SKIN
      allocate( rxyz_old(nat_max,3) )
      allocate( nlist(1:nat_max),list(1:nat_max,nabormx) )
      allocate( cnlist(1:nat_max),clist(1:nat_max,nabormx) )
   END SUBROUTINE

   SUBROUTINE NEW_VLIST(ifirst,ilast)
!     makes the Verlet list using the neighbour list
      USE nlist_mod
!!USE sort_mod
      integer,intent(in):: ifirst,ilast
      real(wp):: r(NDIM),r2
      integer,parameter:: nl_max = 2
      integer:: i,j,ic,jj,nc,nl,nn,neighbors((2*nl_max+1)**3)
      nl = nlayers_ll(rv - tiny(1.0))
      print*,'nl ===============',nl 
      if(nl > nl_max) STOP ' NEW_VLIST: nl > nl_max'
      nlist(ifirst:ilast)=0
      rxyz_old(ifirst:ilast,1:NDIM) = rxyz(ifirst:ilast,1:NDIM)
      call NEW_NLIST(ifirst,ilast)
      print*,'vlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlistvlist'
      do i = ifirst,ilast
      ic = CELL(rxyz(i,1:NDIM))  ! link list 
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
         nc = neighbors(jj)
         if (nc == 0)cycle cell_loop    
         j = hoc(nc)
         do while(j /= 0)
            if (j /= i .AND. j > i) then
            r(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
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
      if (ANY(NLIST(ifirst:ilast) >= nabormx) ) STOP 'NLIST(:) >= nabormx'
!!do i = ifirst,ilast
!!   call qsort(nlist(i),list(i,1:nlist(i)),cclst)
!!   list(i,1:nlist(i))=cclst(1:nlist(i))
!!end do
   END SUBROUTINE NEW_VLIST

   logical function update_vl(i)
      integer,intent(in):: i
      real(wp):: r(3)
      r(1:NDIM) = rxyz(i,1:NDIM)-rxyz_old(i,1:NDIM)
      call pbc(r)
      update_vl = (dot_product(r,r) > SKIN2*0.25_wp)
   end function

   FUNCTION UPDATE2_vl(ifirst,ilast)
!     decides whether the list needs to be reconstructed
      logical:: UPDATE2_vl
      integer,intent(in):: ifirst,ilast
      real(wp):: dispmx
!     a conservative test of the list skin crossing
      dispmx = MAXVAL(ABS( rxyz(ifirst:ilast,1:NDIM) - rxyz_old(ifirst:ilast,1:NDIM) ))
      dispmx = 2.0_wp*SQRT( 3.0_wp*dispmx**2 )
      UPDATE2_vl = ( dispmx > skin )
   END FUNCTION

END MODULE verlet_list_mod


!!>include 'find_atoms_mod_TEOS.f90'

MODULE find_atoms_mod
   USE precision_mod
   USE list_mod

   implicit none
   public:: find_atoms, list_atoms_in_sphere

CONTAINS
   SUBROUTINE find_atoms(iat,ifirst,ilast,atomtype,rad,nat,atomlst)
      use nlist_mod
      use atom_types_mod
      use coordinates_mod
      integer,intent(in):: iat,ifirst,ilast,atomtype
      real(wp),intent(in):: rad
      integer,intent(out):: nat,atomlst(:)
      integer:: ic,nc,j,jj,nl
      real(wp):: dr(NDIM),r2
      integer,parameter:: nl_max = 3
      integer:: neighbors((2*nl_max+1)**3),nn

! Look for atoms of type atomtype within a distance rad of ia
! and atom number between ifirst & ilast
      nl = nlayers_ll(rad - tiny(1.0))
      if (nl > nl_max) STOP 'find_atoms: nl > nl_max'
      r2 = rad*rad

      nat = 0
      atomlst = 0
      ic = CELL(rxyz(iat,1:NDIM))
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
      nc = neighbors(jj)
      if (nc == 0) cycle cell_loop
      j = HOC(nc)
      do while (j /= 0)
         if (j == iat) GOTO 100
         if (atom(j) /= atomtype) GOTO 100
         if (j < ifirst .or. j > ilast) GOTO 100
         dr = rxyz(j,1:NDIM) - rxyz(iat,1:NDIM)
         call pbc(dr)
         if (dot_product(dr,dr) <= r2) then
            nat = nat + 1
            atomlst(nat) = j
         end if
100      j = LL(j)
      end do
      end do cell_loop
   END SUBROUTINE find_atoms


   SUBROUTINE list_atoms_in_sphere(r0,rad,nat,atomlst)
      USE nlist_mod
      USE coordinates_mod
      real(wp),intent(in):: r0(NDIM),rad
      integer,intent(out):: nat,atomlst(:)
      integer:: ic,nc,j,jj,nl
      real(wp):: dr(NDIM),r2
      integer,parameter:: nl_max = 4
      integer:: neighbors((2*nl_max+1)**3),nn
! Look for atoms withinin a distance rad of point r0
      nl = nlayers_ll(rad - tiny(1.0))
      if (nl > nl_max) STOP 'list_atoms_in_sphere: nl > nl_max'
      r2 = rad*rad
      nat = 0
      atomlst = 0
      nl = nlayers_ll(rad)
      ic = CELL(r0)
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
      nc = neighbors(jj)
      if (nc == 0) cycle cell_loop
      j = HOC(nc)
      do while (j /= 0)
         dr = rxyz(j,1:NDIM) - r0
         call pbc(dr)
         if (dot_product(dr,dr) <= r2) then
            nat = nat + 1
            atomlst(nat) = j
         end if
100      j = LL(j)
      end do
      end do cell_loop
   END SUBROUTINE list_atoms_in_sphere

   SUBROUTINE find_glassbond_oxygen(iSi,iOgx)
     USE global_vars_mod
     USE connectivity_mod
     USE rand_mod
     USE list_mod
     integer,intent(in)::iSi
     integer::iSi1,iSi2,i,j,k,iO
     integer,intent(out)::iOgx
     type(mylist)::sio_lst
     real(wp),parameter:: rad = 4.0_wp/angstrom
 iSi1 = iSi
!print*,'insidefindglassbondoxygen'
! print*,'iSi1 =',iSi1, proximity(iSi1,:),atom_name(atomfn(proximity(iSi1,:)))

       call new_list(sio_lst,100)
       do j = 1, ncmax(atom(iSi1))
          k = proximity (iSi1,j)
          if (k == 0) cycle
          if (atom(k) == iOxygen) then
             iSi2 = proximity(k,1)
             if (iSi2 == iSi1) iSi2 = proximity(k,2)
             exit
          else
             cycle
          end if
       end do
       call find_atoms(iSi1,nfixed+1,natom,ioxygen,rad,sio_lst%n,sio_lst%ind)
              i = 1
         do
            if (i > sio_lst%n) exit
            iO = sio_lst%ind(i)
            if ((proximity(iO,1) == iSi1).or.(proximity(iO,2) == iSi1).or.&
                (proximity(iO,1) == iSi2).or.(proximity(iO,2) == iSi2)) then
               sio_lst%ind(i) = sio_lst%ind(sio_lst%n)
               sio_lst%n = sio_lst%n - 1
            else
               i = i + 1
            end if
         end do
         call rand_from_list(SiO_lst,iO)
         iOgx = iO
         call delete_list(SiO_lst)

END SUBROUTINE find_glassbond_oxygen

END MODULE find_atoms_mod

!!>include 'molecule_type_mod.f90'

MODULE molecule_type_mod
   USE precision_mod
   USE atom_types_mod
   USE teos_atomlst_mod
   USE connectivity_mod
   implicit none
!
CONTAINS

SUBROUTINE teos_type(iSi,moltype)
   USE find_atoms_mod
   integer,intent(in):: iSi
   integer,intent(out):: moltype
   integer:: nOc,nOHt,nOx,nfree,nOxfree
   integer:: i,j,iOx,iOgx
!print*, ' == === === === === === === === = == teos_type SUBROUTINE '
!
      nOx = 0
      nOc = 0
      nOHt = 0
      nfree = 0
      nOxfree = 0
      do i = 1, ncmax(atom(iSi))
         j = proximity(iSi,i)
         if (j == 0) then
            nfree = nfree + 1
            cycle
         end if
            select case (atom(j))
            case(iOxygen)
            nOx = nOx + 1
            iOx = j
            case(iOxygenC)
            nOc = nOc + 1
            case(iOxygenH)
            nOHt = nOHt + 1
            case(iOxygenF)
            nOxfree = nOxfree + 1
            end select
      end do
         if (nOc == 3 .and. nOx == 1) then
             moltype = iSiOC3
         else if (nOc == 2 .and. nOHt == 1 .and. nOx == 1) then
            call find_glassbond_oxygen(iSi,iOgx)
            if (iOgx == 0) then
             moltype = iSiOHOC2
            else
             moltype = iSiOHOC2_SiO
            end if
         else if (nOc == 1 .and. nOHt == 2 .and. nOx == 1) then
            call find_glassbond_oxygen(iSi,iOgx)
            if (iOgx == 0) then
             moltype = iSiOH2OC
            else
             moltype = iSiOH2OC_SiO
            end if
         else if (nOc == 1 .and. nOx == 3) then
             moltype = iSiOx3OC1
         else if (nOc == 2 .and. nOx == 2) then
             moltype = iSiOx2OC2
         else if (nOHt == 1 .and. nOx == 3) then
             moltype = iSiOx3OH1
         else if (nOc == 1 .and. nOHt == 1 .and. nOx == 2) then
             moltype = iSiOx2OCOH
         else if (nOx == 4) then
             moltype = iSiOx4
         else if (nOHt == 3.and. nOx == 1) then
             moltype = iOSiOH3
         else if (nOHt == 2.and. nOx == 2) then
             moltype = iSiOx2OH2
         end if

 END SUBROUTINE teos_type
 !
 SUBROUTINE teos_type1(iSi,moltype)
!
    integer,intent(in):: iSi
    integer,intent(out):: moltype
    integer:: nOc,nOHt,nOx,nfree,nOxfree
    integer:: i,j,iOx
!
    nOx = 0
    nOc = 0
    nOHt = 0
    nfree = 0
    nOxfree = 0
    do i = 1, ncmax(atom(iSi))
       j = proximity(iSi,i)
       if (j == 0) then
           nfree = nfree + 1
           cycle
       end if
          select case (atom(j))
          case(iOxygen)
             nOx = nOx + 1
             iOx = j
          case(iOxygenC)
             nOc = nOc + 1
          case(iOxygenH)
             nOHt = nOHt + 1
          case(iOxygenF)
             nOxfree = nOxfree + 1
          end select
    end do
       if (nOc == 3 .and. nOx == 1) then
          moltype = iSiOC3
       else if (nOc == 2 .and. nOHt == 1 .and. nOx == 1) then
          moltype = iSiOHOC2
       else if (nOc == 1 .and. nOHt == 2 .and. nOx == 1) then
          moltype = iSiOH2OC
       else if (nOc == 1 .and. nOx == 3) then
          moltype = iSiOx3OC1
       else if (nOc == 2 .and. nOx == 2) then
          moltype = iSiOx2OC2
       else if (nOHt == 1 .and. nOx == 3) then
          moltype = iSiOx3OH1
       else if (nOc == 1 .and. nOHt == 1 .and. nOx == 2) then
          moltype = iSiOx2OCOH
       else if (nOx == 4) then
          moltype = iSiOx4
       else if (nOHt == 3 .and. nOx == 1) then
          moltype = iOSiOH3
       else if (nOHt == 2 .and. nOx == 2) then
          moltype = iSiOx2OH2
       else if (nOc == 1 .and. nOHt == 3) then
          moltype = iSiOcOH3
       else
          print *,nOHt,nOx,nOxfree,nOc
          call write_atom_info(6,i)
          STOP 'unknown Si'
       end if
   END SUBROUTINE teos_type1

END MODULE molecule_type_mod


!!>include 'charges2_mod.f90'
MODULE charges_mod
      USE precision_mod
      USE atom_types_mod
      USE teos_atomlst_mod
      USE constants_mod
      USE molecule_type_mod
      implicit none
      real(wp),allocatable,public:: charge(:)
      real(wp),allocatable,public:: Ox_gas_charge(:,:),N2_gas_charge(:,:),CO2_gas_charge(:,:),GAS_charge(:,:),h2o_charge(:,:)
      PUBLIC:: assign_charge, qsi
!
CONTAINS

!   PURE FUNCTION qsi(n)
!      real(wp):: qsi
!      integer,intent(in):: n
!      qsi = (n-4)*qi(iOxygen)*0.5_wp - n*(qi(iOxygenH) + qi(iHydrogen))
!   END FUNCTION qsi

   FUNCTION qsi(moltype)
      real(wp):: qsi
      integer,intent(in):: moltype
      real(wp):: q_OC2H5,qO,qOh
      !q_OC2H5 =  qi(iOxygenC) + qi(iCarbon2) + qi(iCarbon3) + 2*qi(iHydrogenC2)+3*qi(iHydrogenC3)
      !q_OC2H5 =  qi(iOxygenC) + qi(iCarbon2) + qi(iCarbon3) + 5*(qi(iHydrogenC))
      q_OC2H5 =  qi(iOxygenC) + qi(iCarbonH2) + qi(iCarbonH3)
      qOh = qi(iOxygenH)+qi(iHydrogen)
      qO = qi(iOxygen)
         select case (moltype)
         case(iSiOC3)
            qsi = -qO*0.5_wp - 3*q_OC2H5
         case(iSiOHOC2)
            qsi = -qO*0.5_wp - qOh - 2*q_OC2H5
         case(iSiOH2OC)
            qsi = -qO*0.5_wp - 2*q_OC2H5
         case(iSiOx3OC1)
            qsi = -qO*1.5_wp - q_OC2H5
         case(iSiOx2OC2)
            qsi = -qO - 2*q_OC2H5
         case(iSiOx3OH1)
            qsi = -qO*1.5_wp - qOh
         case(iSiOx2OCOH)
            qsi = -qO - qOh - q_OC2H5
         case(iOSiOH3)
            qsi = -qO*0.5_wp - 3*qOh
         case(iSiOx4)
            qsi = -qO*2
         case(iSiOx2OH2)
            qsi = -qO - 2*qOh
         case(iSiOcOH3)
            qsi = - 3*qOh - q_OC2H5
         case default
            !print *,'moltype = ',moltype
            stop 'unknown moltype'
         end select
   END FUNCTION qsi


   SUBROUTINE assign_charge(ifirst,ilast)
      USE connectivity_mod
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      real(wp):: qO,qoh,qoc,qc2,qc3,qhc,qch2,qch3
      integer:: i,mtyp
      qoh = qi(iOxygenH)+qi(iHydrogen)
      qO = qi(iOxygen)
      qoc = qi(iOxygenC)
      qc2 = qi(iCarbon2)
      qc3 = qi(iCarbon3)
      qhc = qi(iHydrogenC)
      qch2 = qi(iCarbonH2)
      qch3 = qi(iCarbonH3)
      do i = ifirst,ilast
         select case(atom(i))
         case(iSilicon)
            call teos_type1(i,mtyp)
            charge(i) = qsi(mtyp)
            !print *,i,atom(i),atom_name(atom(i))
         case(iOxygen)
            charge(i) = qO
         case(iOxygenH)
            charge(i) = qoh
         case(iOxygenC)
            charge(i) = qoc
         case(iCarbon2)
            charge(i) = qc2
         case(iCarbon3)
            charge(i) = qc3
         case(iCarbonH2)
            charge(i) = qch2
         case(iCarbonH3)
            charge(i) = qch3
         case(iHydrogenC)
            STOP 'no explicit Hc used'
         case(iHydrogen)
            STOP 'no explicit H used'
            !charge(i) = qh
         case default
            print *,i,atom(i),atom_name(atom(i))
            stop 'assign_charge: unknown atom type'
         end select
      end do
   END SUBROUTINE assign_charge

END MODULE charges_mod



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

   SUBROUTINE repul_force_List(nl,list,en)
      USE coordinates_mod
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      integer,intent(in):: nl,list(:)
      real(wp),intent(out):: en
      real(wp):: rr(NDIM),r2,r(NDIM)
      integer:: M,L,iL,ic,nc,inn,nn,neigbors(27)
!
      en = 0.0_wp
      list_loop: do iL = 1, nl
      L = list(il)
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
      end do list_loop
      en = en*0.5_wp*ff
   END SUBROUTINE repul_force_List


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
            if (r2 <= rden2) then
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
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            if (r2 <= rden2) then
               en = en + (r2 - rden2)**2
            end if
         end do atom_loop
      end do outer_atom_loop
      repul_energy_VL = en*0.5_wp*ff
   END FUNCTION repul_energy_VL

END MODULE repul_energy_mod

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
   END SUBROUTINE ran_point_in_sphere

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
   END SUBROUTINE ran3sph

   SUBROUTINE ran3sph2(x,y,z)
      real(wp),parameter::twopi = 6.283185307179586476925286766559005768394_wp
      real(wp),intent(out)::x,y,z
      real(wp)::s,tmp
      z = 2.0_wp*rand() - 1.0_wp
      tmp = rand()*twopi
      s = sqrt(1.0_wp - z*z)
      x = s*cos(tmp)
      y = s*sin(tmp)
   END SUBROUTINE ran3sph2

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
   END SUBROUTINE ran_point_on_sphere

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
   END SUBROUTINE ran4sph

END MODULE ran_point_sphere


!!>include 'insert_mod_2.f90'

MODULE insert_mod
    USE precision_mod, only: wp
    USE global_vars_mod
    USE nlist_mod
    USE atom_types_mod
    USE coordinates_mod
    USE constants_mod
    USE rand_mod
    USE quat2mat_mod
    USE ran_point_sphere
    implicit none

CONTAINS

   SUBROUTINE insert_mol_GAS(n_mol,n_atom,mol_vec,real_xyz,ntrial,rad1,rlower,rupper,fover,self_overlap_check)
      integer,intent(in):: ntrial
      integer,intent(inout):: n_mol,n_atom
      real(wp),intent(in):: mol_vec(:,:),rad1(n_atom)
      real(wp),intent(in):: rlower(:),rupper(:)
      real(wp),intent(out):: real_xyz(:,:,:)
      real(wp),intent(in),optional:: fover
      logical,intent(in),optional:: self_overlap_check
      real(wp):: aa(3,3),q(4),x,y,z,fover0
      real(wp):: mol_vec1(size(mol_vec,1),size(mol_vec,2))
      integer:: it,i,im,j,n_mol_requested
      logical:: self_overlap_check0
!
      if (present(fover)) then
         fover0 = fover
      else
         fover0 = 1.0_wp
      end if
!
      if (present(self_overlap_check)) then
         self_overlap_check0 = self_overlap_check
      else
         self_overlap_check0 = .FALSE.
      end if
!
      n_mol_requested = n_mol
      n_mol = 0
      mol_loop: do im = 1, n_mol_requested
      insert_loop: do it = 1,ntrial

         x = rlower(1) + rand()*(rupper(1) - rlower(1))
         y = rlower(2) + rand()*(rupper(2) - rlower(2))
         z = rlower(3) + rand()*(rupper(3) - rlower(3))

         !call random_rotation(aa)
         call ran4sph(q)
         call quat2mat(q,aa)
         do i = 1,n_atom
            mol_vec1(:,i) = matmul(aa,mol_vec(:,i))
         end do

         mol_vec1(1,1:n_atom) = mol_vec1(1,1:n_atom) + x
         mol_vec1(2,1:n_atom) = mol_vec1(2,1:n_atom) + y
         mol_vec1(3,1:n_atom) = mol_vec1(3,1:n_atom) + z
         do i = 1,n_atom
            call pbc(mol_vec1(1:NDIM,i))
         end do

! first check if there is a large overlap
         if (overlap(fover0,rad1,mol_vec1,n_atom)) cycle insert_loop
         if (self_overlap_check0) then
            if (overlap2(fover0,rad1,mol_vec1,n_atom,n_mol*n_atom,&
                         ratm=reshape(real_xyz(1:n_mol,1:n_atom,:),(/n_mol*n_atom,3/)),&
                         radatm=reshape(spread(rad1,dim=1,ncopies=n_mol),(/n_mol*n_atom/)) )) cycle insert_loop
         end if

         n_mol = n_mol + 1
         do j = 1,n_atom
            real_xyz(n_mol,j,:) = mol_vec1(:,j)
         end do
         write(*,*) im,n_mol
         cycle mol_loop

      end do insert_loop
      end do mol_loop
      if (n_mol == 0) n_atom = 0
   END SUBROUTINE insert_mol_GAS


   FUNCTION overlap(fr,rad,mol_vec1,n_atom)
      logical:: overlap
      integer,intent(in):: n_atom
      integer:: k,jj,j,ic,nc
      real(wp):: dr(NDIM),rad(:),fr,mol_vec1(:,:),sigj2
      integer:: neighbors(125),nn
      overlap = .false.
      probe_atom_loop: do k = 1,n_atom
         ic = CELL( mol_vec1(1:NDIM,k) )
         call NEIGCELL(ic,2,nn,neighbors)
         cell_loop: do jj = 1,nn
            nc = neighbors(jj)
            if (nc == 0) cycle cell_loop
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               sigj2 = sigLJ_2(atom(j))
               if (atom(j) == iSilicon) then
                  sigj2 = 1.0_wp/Angstrom
               end if
               dr(1:NDIM) = mol_vec1(1:NDIM,k) - rxyz(j,1:NDIM)
               call pbc(dr)
               if (dot_product(dr,dr) < (fr*(rad(k) + sigj2))**2) then
                  overlap = .true.
                  RETURN
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
      end do probe_atom_loop
   END FUNCTION overlap


   FUNCTION overlap2(fr,rad,mol_vec1,n_atom,nat,ratm,radatm)
      logical:: overlap2
      integer,intent(in):: n_atom,nat
      real(wp),intent(in):: fr,rad(:),mol_vec1(:,:),ratm(:,:),radatm(:)
      integer:: k,j
      real(wp):: dr(NDIM)
      overlap2 = .false.
      probe_atom_loop: do k = 1,n_atom
         do j = 1,nat
            dr(1:NDIM) = mol_vec1(1:NDIM,k) - ratm(j,1:NDIM)
            call pbc(dr)
            if (dot_product(dr,dr) < (fr*(rad(k) + radatm(j)))**2) then
               overlap2 = .true.
               RETURN
            end if
         end do
      end do probe_atom_loop
   END FUNCTION overlap2

END MODULE insert_mod

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

   SUBROUTINE rotate_inverseZaxis_mol(rxyz,r1,n_by_3)
      real(wp),intent(inout):: rxyz(:,:)
      real(wp),intent(in):: r1(:)
      logical,intent(in):: n_by_3
      real(wp):: at(3,3),aa(3,3)
      integer:: i,n1,n2
      n1 = size(rxyz,1)
      n2 = size(rxyz,2)
      call zaxis2vect(r1,aa)
      at = transpose(aa)
      if (n_by_3) then
         do i = 1,n1
            rxyz(i,:) = matmul(at,rxyz(i,:))
         end do
         rxyz(1:n1,2:3) = -rxyz(1:n1,2:3)
      else
         rxyz(:,1:n2) = matmul(at,rxyz(:,1:n2))
         rxyz(2:3,1:n2) = -rxyz(2:3,1:n2)
      end if
   END SUBROUTINE rotate_inverseZaxis_mol

END MODULE rotate_axis_mod

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
   END SUBROUTINE tetra_coords

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
   END SUBROUTINE tetra_coords5

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
   END SUBROUTINE SiOH4_coords

END MODULE tetra_coords_mod


!!>include 'probe_mol_mod.f90'

MODULE probe_mol_mod
   USE precision_mod
   USE constants_mod
   USE atom_types_mod
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
      integer,allocatable:: nc(:),conect(:,:)
   END TYPE probe_mol
!

CONTAINS

   SUBROUTINE new_probe_mol(P,n,radius,atomtype,flexible)
      type(probe_mol),intent(out):: P
      integer,intent(in):: n
      real(wp),intent(in),optional:: radius
      integer,intent(in),optional:: atomtype
      logical,intent(in),optional:: flexible
      P%n = n
      allocate(P%r(3,n),P%atom(n),P%rad(n),P%q(n))
      if (present(flexible)) then
         if (flexible) then
            allocate(P%nc(n),P%conect(n,4))
            P%nc = 0
            P%conect = 0
         end if
      end if
      P%r = 0.0_wp
      P%atom = 0
      P%rad = 0.0_wp
      if (present(radius)) P%rad = radius
      if (present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE


   SUBROUTINE delete_probe_mol(P)
      type(probe_mol),intent(inout):: P
      P%n = 0
      if (allocated(P%r)) deallocate(P%r)
      if (allocated(P%atom)) deallocate(P%atom)
      if (allocated(P%rad)) deallocate(P%rad)
      if (allocated(P%q)) deallocate(P%q)
      if (allocated(P%nc)) deallocate(P%nc)
      if (allocated(P%conect)) deallocate(P%conect)
   END SUBROUTINE delete_probe_mol


   SUBROUTINE WRITE_PROBE_MOL(iu,P,remark)
      USE atom_types_mod
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      character(*),intent(in):: remark
      integer:: i,j
      write(iu,*) P%n
      write(iu,*) trim(remark)
      do i = 1,P%n
!        write(iu,'(a6,i5,a4,2x,a3,i6,4x,3f8.3')'HETATM',i,atom_name(P%atom(i)), &
!                                               '   ',i,P%r(1:3,i)*angstrom
         write(iu,'(a3,5(1x,f11.5))') atom_name(P%atom(i)),P%r(1:3,i)*Angstrom,P%rad(i),P%q(i)
      end do
      if (allocated(P%conect)) then
         do i = 1,P%n
            write(iu,'(2(i4,4x),4i4)') i,P%nc(i),(P%conect(i,j),j=1,P%nc(i))
         end do
      end if
      write(iu,*)
   END SUBROUTINE

END MODULE probe_mol_mod

!!>include 'probe_mol_init_mod_TEOS.f90'

MODULE probe_mol_init_mod
   USE precision_mod
   USE probe_mol_mod
!
   type(probe_mol),save:: probe_SiOH4, probe_SiO4, probe_tip3p,probe_SiO4Q,probe_SiO3
   type(probe_mol),save:: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet,probe_TEOS
!
CONTAINS

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types_mod
      USE constants_mod
      USE charges_mod
      USE rotate_axis_mod
!      USE coordinates_mod, only: read_pdb_mol
      real(wp):: r1(3)
      integer:: i
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigLJ_H2O*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom = (/ iOw, iHw, iHw /)
      probe_TIP3P%rad = sigLJ_2(probe_TIP3P%atom)
      probe_TIP3P%q = (/q_O_H2O, -0.5_wp*q_O_H2O, -0.5_wp*q_O_H2O /)
  
! Si(OH)4
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO,bondl_OH,angle_SiOH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigLJ_2(probe_SiOH4%atom)
      ! The Silicon atom should have a radius >= 1.0 A
      probe_SiOH4%rad(1) = max(1.0_wp/Angstrom,probe_SiOH4%rad(1))
      probe_SiOH4%q = qi(probe_SiOH4%atom)
      probe_SiOH4%q(1) = -4.0_wp*q_OH_sil

! SiO4, just the 4 OH Oxygens
      call new_probe_mol(probe_SiO4,4,sigLJ_OH*0.5_wp)
      call tetra_coords(bondl_SiO,probe_SiO4%r)
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigLJ_2(probe_SiO4%atom)

! SiO4Q, Si & 4 united atom OH Oxygens with charges & connectivity
      call new_probe_mol(probe_SiO4Q,5)
      call tetra_coords5(bondl_SiO,probe_SiO4Q%r)
      probe_SiO4Q%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4Q%rad = sigLJ_2(probe_SiO4Q%atom)
      ! The Silicon atom should have a radius >= 1.0 A
      probe_SiOH4%rad(1) = max(1.0_wp/Angstrom,probe_SiOH4%rad(1))
      probe_SiO4Q%q(2:5) = q_OH_sil
      probe_SiO4Q%q(1) = -4.0_wp*q_OH_sil
! Debugging: check molecule is electroneutral
!if(sum(probe_SiO4Q%q)/=0.0_wp)then
!print *,probe_SiO4Q%q
!print *,sum(probe_SiO4Q%q)
!stop 'probe_SiO4Q%q/=0'
!end if
      allocate(probe_SiO4Q%nc(5),probe_SiO4Q%conect(5,4))
      probe_SiO4Q%nc = (/ 4,1,1,1,1/)
      probe_SiO4Q%conect(1,:) = (/ 2,3,4,5 /)
      probe_SiO4Q%conect(2,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%conect(3,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%conect(4,:) = (/ 1,0,0,0 /)
      probe_SiO4Q%conect(5,:) = (/ 1,0,0,0 /)
!
! SiO3, Si & 3 OH Oxygens
      call new_probe_mol(probe_SiO3,4)
      probe_SiO3%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO3%rad = sigLJ_2(probe_SiO3%atom)
      allocate(probe_SiO3%nc(4),probe_SiO3%conect(4,4))
      probe_SiO3%nc = (/ 3,1,1,1 /)
      probe_SiO3%conect(1,:) = (/ 2,3,4,0 /)
      probe_SiO3%conect(2,:) = (/ 1,0,0,0 /)
      probe_SiO3%conect(3,:) = (/ 1,0,0,0 /)
      probe_SiO3%conect(4,:) = (/ 1,0,0,0 /)
      probe_SiO3%r(:,1) = probe_SiO4Q%r(:,1)
      probe_SiO3%r(:,2) = probe_SiO4Q%r(:,3)
      probe_SiO3%r(:,3) = probe_SiO4Q%r(:,4)
      probe_SiO3%r(:,4) = probe_SiO4Q%r(:,5)

! Si(OC2H5)4
 !     call new_probe_mol(probe_TEOS,33,flexible=.true.)
      call new_probe_mol(probe_TEOS,13)
      open(771,file='teos2.ent',status='old')
      call read_pdb_mol1(771,probe_TEOS%N,probe_TEOS%r,probe_TEOS%atom)
      if (probe_TEOS%atom(1) /= iSilicon) STOP 'incorrect TEOS molecule'

      forall(i = 1:probe_TEOS%N) probe_TEOS%r(:,i) = probe_TEOS%r(:,i) - probe_TEOS%r(:,1)
      r1 = probe_TEOS%r(:,11) - probe_TEOS%r(:,1)
      r1 = r1/sqrt(dot_product(r1,r1))
      call rotate_inverseZaxis_mol(probe_TEOS%r,r1,n_by_3 = .false.)
print*,'========================================================'
write(*,*) probe_TEOS%N
write(*,*)'TEOS'
do i = 1,probe_TEOS%N
write(*,'(a,3f12.6)') atom_name(probe_TEOS%atom(i)),probe_TEOS%r(:,i)
end do
print*,'========================================================'
      probe_TEOS%rad = sigLJ_2(probe_TEOS%atom)
      ! The Silicon atom should have a radius >= 1.0 A
      probe_TEOS%rad(1) = max(1.0_wp/Angstrom,probe_TEOS%rad(1))

      probe_TEOS%q = qi(probe_TEOS%atom)
print*,'probe_TEOS%q =',probe_TEOS%q
      probe_TEOS%q(1) = -sum(probe_TEOS%q(2:probe_TEOS%n))
print*,'probe_TEOS%q(1) =',probe_TEOS%q(1)
print *,sum(probe_TEOS%q)
      allocate(probe_TEOS%nc(13),probe_TEOS%conect(13,4))
      probe_TEOS%nc = (/ 4,2,2,1,2,2,1,&
                           2,2,1,2,2,1 /)

  call  connect_TEOS2(probe_TEOS%conect,probe_TEOS%atom)
  print*, 'probe_TEOS%conect',probe_TEOS%conect(:,:)

! O2
      call new_probe_mol(probe_O2,2)
      probe_O2%atom = (/ iO_O2, iO_O2 /)
      probe_O2%rad = sigLJ_2(probe_O2%atom)
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_O2*0.5_wp /)
      probe_O2%q = qi(probe_O2%atom)
! N2
      call new_probe_mol(probe_N2,2)      
      probe_N2%atom = (/ iN_N2, iN_N2 /)
      probe_N2%rad = sigLJ_2(probe_N2%atom)
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp,-bondl_N2*0.5_wp /)
      probe_N2%q = qi(probe_N2%atom)
! CO2
      call new_probe_mol(probe_CO2,3)
!     probe_CO2%rad(1) = sig_C_CO2*0.5_wp     
      probe_CO2%atom = (/ iO_CO2, iC_CO2, iO_CO2 /)
      probe_CO2%rad = sigLJ_2(probe_CO2%atom)
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp,-bondl_CO2 /)
      probe_CO2%q = qi(probe_CO2%atom)
!
! Write out details of the probe molecules
! Si(OH)4
      call write_probe_mol(6,probe_SiOH4,'SiOH4')
      print '(a,f12.6/)','Total SiOH4 charge = ',sum(probe_SiOH4%q)
      call write_probe_mol(6,probe_TIP3P,'H2O TIP3P')
      call write_probe_mol(6,probe_SiO4,'SiO4 ... just OH')
      call write_probe_mol(6,probe_SiO4Q,'SiO4')
      call write_probe_mol(6,probe_SiO3,'SiO3 ... just OH')
      call write_probe_mol(6,probe_O2,'O2')
      call write_probe_mol(6,probe_N2,'N2')
      call write_probe_mol(6,probe_CO2,'CO2')
      call write_probe_mol(6,probe_TEOS,'TEOS = Si(OC2H5)4')
      print '(a,f12.8/)','Total TEOS charge = ',sum(probe_TEOS%q)
      print '(a,f12.8/)','Total CO2 charge = ',sum(probe_CO2%q)
      print '(a,f12.8/)','Total O2 charge = ',sum(probe_O2%q)
      print '(a,f12.8/)','Total N2 charge = ',sum(probe_N2%q)
!
   END SUBROUTINE init_probe_mols


   SUBROUTINE h2o_coords(bondl,ang,r)
      real(wp),intent(in):: bondl,ang
      real(wp),intent(out):: r(3,3)
! Coordinates for water molecule
      r(1:3,1) = 0.0_wp
      r(1:3,2) = (/ 0.0_wp, 0.0_wp, bondl /)
      r(1:3,3) = (/ 0.0_wp, bondl*sin(ang), bondl*cos(ang) /)
   END SUBROUTINE
   SUBROUTINE read_pdb_mol1(iu,N,rxyz,atom)
      integer,intent(in):: iu,N
      real(wp),intent(out):: rxyz(:,:)
      integer,intent(out):: atom(:)
      integer:: i,itmp
      character(6):: ctmp1,ctmp2
      do i = 1,N
         read(iu,*) ctmp1,itmp,ctmp2,itmp,rxyz(1:3,i)
!print*,ctmp1,itmp,ctmp2,itmp,rxyz(1:3,i)
         rxyz(1:3,i) = rxyz(1:3,i)/Angstrom
         atom(i) = name2atom(trim(ctmp2))
!print*, trim(ctmp2)
!print*,atom(i)
      end do
!110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3)
   END SUBROUTINE read_pdb_mol1


   SUBROUTINE connect_TEOS(conect,atype)
      USE atom_types_mod
      integer,intent(out):: conect(:,:),atype(:)
!
!                             5  7
!        H H               2--3--4--8
!      O-C-C-H            /   6  9
!     /  H H             /
!    /   H H            /     13  15
!  Si--O-C-C-H         1--10--11--12--16
!  | \   H H           |\     14  17
!  |  \  H H           | \
!  |   O-C-C-H         |  \    21  23
!  |     H H           |   18--19--20--24
!  |                   |       22  25
!  \  H H              |
!   O-C-C-H            \    29  31
!     H H               26--27--28--32
!                           30  33
!
      conect = 0
      conect(1,1:4) = (/ 2,10,18,26 /)
      conect(2,1:4) = (/ 1,3,0,0 /)
      conect(3,1:4) = (/ 2,4,5,6 /)
      conect(4,1:4) = (/ 3,7,8,9 /)
      conect(5,1:4) = (/ 3,0,0,0 /)
      conect(6,1:4) = (/ 3,0,0,0 /)
      conect(7,1:4) = (/ 4,0,0,0 /)
      conect(8,1:4) = (/ 4,0,0,0 /)
      conect(9,1:4) = (/ 4,0,0,0 /)
      where (conect(2:9,:) /= 0)
         conect(10:17,:) = conect(2:9,:) + 8
         conect(18:25,:) = conect(2:9,:) + 16
         conect(26:33,:) = conect(2:9,:) + 24
      end where
      conect(2,1) = 1
      conect(10,1) = 1
      conect(18,1) = 1
      conect(26,1) = 1
      atype(1) = iSilicon
      atype(2:9) = (/ ioxygenC,icarbon2,icarbon3,iHydrogenC,iHydrogenC,iHydrogenC,iHydrogenC,iHydrogenC /)
      atype(10:17) = atype(2:9)
      atype(18:25) = atype(2:9)
      atype(26:33) = atype(2:9)
   END SUBROUTINE connect_TEOS

   
   SUBROUTINE connect_TEOS2(conect,atype)
      USE atom_types_mod
      integer,intent(out):: conect(:,:),atype(:)
!
!                               
!                           2--3--4
!      O-Ch2-Ch3           /   
!     /                  /
!    /                  /       
!  Si--O-Ch2-Ch3       1--5--6--7
!  | \                 |\     
!  |  \                | \
!  |   O-Ch2-Ch3       |  \    
!  |                   |   8--9--10
!  |                   |       
!  \                   |
!   O-Ch2-Ch3          \    
!                       11--12--13
!                           
!
      conect = 0
      conect(1,1:4) = (/ 2,5,8,11 /)
      conect(2,1:4) = (/ 1,3,0,0 /)
      conect(3,1:4) = (/ 2,4,0,0 /)
      conect(4,1:4) = (/ 3,0,0,0 /)

      where (conect(2:4,:) /= 0)
         conect(5:7,:) = conect(2:4,:) + 3
         conect(8:10,:) = conect(2:4,:) + 6
         conect(11:13,:) = conect(2:4,:) + 9
      end where
      conect(2,1) = 1
      conect(5,1) = 1
      conect(8,1) = 1
      conect(11,1) = 1
      atype(1) = iSilicon
      atype(2:4) = (/ iOxygenC,iCarbonH2,iCarbonH3 /)
      atype(5:7) = atype(2:4)
      atype(8:10) = atype(2:4)
      atype(11:13) = atype(2:4)
   END SUBROUTINE connect_TEOS2

END MODULE probe_mol_init_mod


!!>include 'lj_el_mod_O2_N2_CO2.f90'

MODULE Lj_el_mod
    USE precision_mod
    USE global_vars_mod
    USE coordinates_mod
    USE charges_mod
    USE atom_types_mod
    USE constants_mod
    USE nlist_mod
    implicit none
    save
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: vdwcut(:,:),vdwsf(:,:),vdwcut_mod(:,:),vdwsf_mod(:,:)
    real(wp),parameter:: rcut = 12.0_wp/Angstrom
    real(wp),parameter:: rcut2 = rcut**2
    real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
    real(wp),parameter:: rcut7 = rcut6/rcut
    real(wp),parameter:: rcut13 = rcut7*rcut7*rcut
!   real(wp),parameter:: rcut_1 = 3.32_wp/Angstrom
!   real(wp),parameter:: rcut_1_2 = rcut_1**2
!   real(wp),parameter:: rcut_1_6 = 1.0_wp/(rcut_1_2*rcut_1_2*rcut_1_2)
!   real(wp),parameter:: rcut_1_7 = 1.0_wp/(((rcut_1**3)**2)*rcut_1)
!   real(wp),parameter:: rcut_1_13 = rcut_1_7*rcut_1_7*rcut_1

CONTAINS

   SUBROUTINE LJ_INIT
      real(wp):: eij,rstij
      integer:: i,j
      allocate( aij(ntyplj,ntyplj) )
      allocate( bij(ntyplj,ntyplj) )
      write(*,*) 'ntyplj = ',ntyplj
      allocate( vdwcut(ntyplj,ntyplj) )
      allocate( vdwsf(ntyplj,ntyplj) )
      allocate( vdwcut_mod(ntyplj,ntyplj) )
      allocate( vdwsf_mod(ntyplj,ntyplj) )
!
      do j = 1,ntyplj    ! lj aij,bij parameters & shifted terms
      do i = 1,ntyplj
         eij = sqrt(epsLJ(i)*epsLJ(j))
         rstij = (sigLJ(i) + sigLJ(j))*0.5_wp
         aij(i,j) = 4.0_wp*(rstij**12)*eij
         bij(i,j) = 4.0_wp*(rstij**6)*eij
         vdwcut(i,j) = rcut6*(aij(i,j)*rcut6 - bij(i,j))
         vdwsf(i,j) = rcut6*(-12.0_wp*aij(i,j)*rcut6 + 6.0_wp*bij(i,j))/rcut
!        vdwcut_mod(i,j) = rcut_1_6*(aij(i,j)*rcut_1_6-bij(i,j))
!        vdwsf_mod(i,j) =  rcut_1_6*(-12.0_wp*aij(i,j)*rcut_1_6+6.0*bij(i,j))/rcut_1
!        write(*,*) i,j,aij(i,j),bij(i,j),vdwcut(i,j),vdwsf(i,j)
      end do
      end do
   END SUBROUTINE LJ_INIT


   SUBROUTINE FORCE_ENERGY_LJ_EL_O2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         image_Ox_xyz(M,1,:) = Ox_xyz(M,1,:)
         rr(1:NDIM) = Ox_xyz(M,2,:) - Ox_xyz(M,1,:)
         call pbc(rr)
         image_Ox_xyz(M,2,:) = Ox_xyz(M,1,:) + rr(:)
         image_Ox_xyz(M,3,:) = (image_Ox_xyz(M,1,:) + image_Ox_xyz(M,2,:))/2.0_wp

       do n = 1,2
         rr(1:NDIM) = Ox_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2*r2*r2)
            rr13 = rr6*rr6/r1
            rr7 = rr6/r1
            dele = rr6*(aij(atom(L),Ox_atom(M,n))*rr6 - bij(atom(L),Ox_atom(M,n))) &
                 - vdwcut(atom(L),Ox_atom(M,n)) - vdwsf(atom(L),Ox_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),Ox_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),Ox_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            Ox_fxyz(M,n,1:NDIM) = Ox_fxyz(M,n,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,Ox_atom(M,n))*rr6-bij(0,Ox_atom(M,n))) &
!                 - vdwcut_mod(0,Ox_atom(M,n)) - vdwsf_mod(0,Ox_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,Ox_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,Ox_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            Ox_fxyz(M,n,1:NDIM) = Ox_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + dele
!         end if

!         rr(1:NDIM) = image_Ox_xyz(M,n,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = dot_product(rr,rr)
!         r1 = sqrt(r2)
!            df = -charge(L)*Ox_gas_charge(M,n)*rr(1:NDIM)/(r2*r1)
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            Ox_fxyz(M,n,1:NDIM) = Ox_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,n)/r1
!!

      end do
!
!         rr(1:NDIM) = image_Ox_xyz(M,3,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = dot_product(rr,rr)
!         r1 = sqrt(r2)
!         df = -charge(L)*Ox_gas_charge(M,3)*rr(1:NDIM)/(r2*r1)
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            Ox_fxyz(M,1,1:NDIM) = Ox_fxyz(M,1,1:NDIM) - 0.5_wp*df
!            Ox_fxyz(M,2,1:NDIM) = Ox_fxyz(M,2,1:NDIM) - 0.5_wp*df
!         ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,3)/r1
!
!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,Ox_atom(M,2))*rr6-bij(0,Ox_atom(M,2))) &
!                 - vdwcut_mod(0,Ox_atom(M,2)) - vdwsf_mod(0,Ox_atom(M,2))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,Ox_atom(M,2))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,Ox_atom(M,2))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            Ox_fxyz(M,1,1:NDIM) = Ox_fxyz(M,1,1:NDIM) - 0.5_wp*df
!            Ox_fxyz(M,2,1:NDIM) = Ox_fxyz(M,2,1:NDIM) - 0.5_wp*df
!            ULJEL = ULJEL + dele
!         end if

      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE FORCE_ENERGY_LJ_EL_O2


   SUBROUTINE FORCE_ENERGY_LJ_EL_N2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         image_N2_xyz(M,1,:) = N2_xyz(M,1,:)
         rr(1:NDIM) = N2_xyz(M,2,:) - N2_xyz(M,1,:)
         call pbc(rr)
         image_N2_xyz(M,2,:) = N2_xyz(M,1,:) + rr(:)
         image_N2_xyz(M,3,:) = (image_N2_xyz(M,1,:) + image_N2_xyz(M,2,:))/2.0_wp

       do n = 1,2
         rr(1:NDIM) = N2_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            dele = rr6*(aij(atom(L),N2_atom(M,n))*rr6 - bij(atom(L),N2_atom(M,n))) &
                 - vdwcut(atom(L),N2_atom(M,n)) - vdwsf(atom(L),N2_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),N2_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),N2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            N2_fxyz(M,n,1:NDIM) = N2_fxyz(M,n,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,N2_atom(M,n))*rr6-bij(0,N2_atom(M,n))) &
!                 - vdwcut_mod(0,N2_atom(M,n)) - vdwsf_mod(0,N2_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,N2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,N2_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            N2_fxyz(M,n,1:NDIM) = N2_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + dele
!         end if
!
!         rr(1:NDIM) = image_N2_xyz(M,n,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = dot_product(rr,rr)
!         r1 = sqrt(r2)
!            df = -charge(L)*N2_gas_charge(M,n)*rr(1:NDIM)/(r2*r1)
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            N2_fxyz(M,n,1:NDIM) = N2_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + charge(L)*N2_gas_charge(M,n)/r1
      end do

!         rr(1:NDIM) = image_N2_xyz(M,3,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = dot_product(rr,rr)
!         r1 = sqrt(r2)
!         df = -charge(L)*N2_gas_charge(M,3)*rr(1:NDIM)/(r2*r1)
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            N2_fxyz(M,1,1:NDIM) = N2_fxyz(M,1,1:NDIM) - 0.5_wp*df
!            N2_fxyz(M,2,1:NDIM) = N2_fxyz(M,2,1:NDIM) - 0.5_wp*df
!         ULJEL = ULJEL + charge(L)*N2_gas_charge(M,3)/r1
!
!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,N2_atom(M,2))*rr6-bij(0,N2_atom(M,2))) &
!                 - vdwcut_mod(0,N2_atom(M,2)) - vdwsf_mod(0,N2_atom(M,2))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,N2_atom(M,2))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,N2_atom(M,2))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            N2_fxyz(M,1,1:NDIM) = N2_fxyz(M,1,1:NDIM) - 0.5_wp*df
!            N2_fxyz(M,2,1:NDIM) = N2_fxyz(M,2,1:NDIM) - 0.5_wp*df
!            ULJEL = ULJEL + dele
!         end if

      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE FORCE_ENERGY_LJ_EL_N2


   SUBROUTINE FORCE_ENERGY_LJ_EL_CO2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         image_CO2_xyz(M,1,:) = CO2_xyz(M,1,:)
         rr(1:NDIM) = CO2_xyz(M,2,:) - CO2_xyz(M,1,:)
         call pbc(rr)
         image_CO2_xyz(M,2,:) = CO2_xyz(M,1,:) + rr(:)
         rr(1:NDIM) = CO2_xyz(M,3,:) - CO2_xyz(M,1,:)
         call pbc(rr)
         image_CO2_xyz(M,3,:) = CO2_xyz(M,1,:) + rr(:)

         do n = 1,3
         rr(1:NDIM) = CO2_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            dele = rr6*(aij(atom(L),CO2_atom(M,n))*rr6 - bij(atom(L),CO2_atom(M,n))) &
                 - vdwcut(atom(L),CO2_atom(M,n)) - vdwsf(atom(L),CO2_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),CO2_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),CO2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            CO2_fxyz(M,n,1:NDIM) = CO2_fxyz(M,n,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,CO2_atom(M,n))*rr6-bij(0,CO2_atom(M,n))) &
!                 - vdwcut_mod(0,CO2_atom(M,n)) - vdwsf_mod(0,CO2_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,CO2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,CO2_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            CO2_fxyz(M,n,1:NDIM) = CO2_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + dele
!         end if
!
!            rr(1:NDIM) = image_CO2_xyz(M,n,:) - rxyz(L,:)
!            call pbc(rr)
!            r2 = dot_product(rr,rr)
!            r1 = sqrt(r2)
!
!            df = -charge(L)*CO2_gas_charge(M,n)*rr(1:NDIM)/(r2*r1)
!            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
!            CO2_fxyz(M,n,1:NDIM) = CO2_fxyz(M,n,1:NDIM) - df
!            ULJEL = ULJEL + charge(L)*CO2_gas_charge(M,n)/r1

         end do

      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE FORCE_ENERGY_LJ_EL_CO2

   SUBROUTINE ENERGY_LJ_EL(pr,sysfirst,syslast,Ulj,Uel)
      USE coordinates_mod
      USE atom_types_mod
      USE charges_mod, only: charge
      USE probe_mol_mod
      type(probe_mol),intent(in):: pr
      integer,intent(in):: sysfirst,syslast
      real(wp),intent(out):: Ulj,Uel
      real(wp):: rr(NDIM),r2,r1,rr6,dele,qL
      integer:: L,M,iaL
      Ulj = 0.0_wp
      Uel = 0.0_wp
!print*,'pr%n=',pr%n
      outer_atom_loop: do L = 1,pr%n
      iaL = pr%atom(L)
!print*,'iaL=',iaL
      qL = pr%q(L)
      atom_loop: do M = sysfirst,syslast
         rr = rxyz(M,1:3) - pr%r(1:3,L)
         call pbc(rr)
         !r2 = dot_product(rr,rr)
         r2 = dot_product(rr,rr)
         r1 = sqrt(r2)
         if (r2 < rcut2 .and. M > 0) then
            rr6 = 1.0_wp/(r2**3)
            dele = rr6*(aij(iaL,atom(M))*rr6-bij(iaL,atom(M))) !&
                 !- uljcut(iaL,atom(M)) - uljsf(iaL,atom(M))*(r1-rcut)
!print*,'dele=',dele
            Ulj = Ulj + dele
!print*,'Ulj=',Ulj
         end if
         Uel = Uel + qL*charge(M)/r1
!if (Uel < -10.0_wp)stop
!print*,'Uel=',Uel
      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE ENERGY_LJ_EL


   SUBROUTINE O2_stretching(en)
!     Force routine for Keating bond stretching
      real(wp),intent(out):: en
      real(wp):: rl1(NDIM),r1sq !,df(NDIM)
      integer:: i
      en = 0.0_wp
      do i = 1,n_O2
        rl1 = Ox_xyz(i,1,:) - Ox_xyz(i,2,:)
        call pbc(rl1)
        r1sq = sqrt(dot_product(rl1,rl1))
        Ox_fxyz(i,1,1:NDIM) = Ox_fxyz(i,1,1:NDIM) - KOO*(r1sq - AOO)*rl1(1:NDIM)/r1sq
        Ox_fxyz(i,2,1:NDIM) = Ox_fxyz(i,2,1:NDIM) + KOO*(r1sq - AOO)*rl1(1:NDIM)/r1sq
        en = en + (r1sq - AOO)**2
      end do
      en = en*0.5_wp*KOO
    END SUBROUTINE O2_stretching

   SUBROUTINE N2_stretching(en)
!     Force routine for Keating bond stretching
      real(wp),intent(out):: en
      real(wp):: r1(NDIM),r1sq !,df(NDIM)
      integer:: i
      en = 0.0_wp
      do i = 1,n_N2
        r1 = N2_xyz(i,1,:) - N2_xyz(i,2,:)
        call pbc(r1)
        r1sq = sqrt(dot_product(r1,r1))
        N2_fxyz(i,1,1:NDIM) = N2_fxyz(i,1,1:NDIM) - KNN*(r1sq - ANN)*r1(1:NDIM)/r1sq
        N2_fxyz(i,2,1:NDIM) = N2_fxyz(i,2,1:NDIM) + KNN*(r1sq - ANN)*r1(1:NDIM)/r1sq
        en = en + (r1sq - ANN)**2
      end do
      en = en*0.5_wp*KNN
   END SUBROUTINE N2_stretching


   SUBROUTINE CO2_stretching(en)
!     Force routine for Keating bond stretching
      real(wp),intent(out):: en(:)
      real(wp):: rl1(NDIM),rl2(NDIM),df1(NDIM),df2(NDIM),r1sq,r2sq,dist
      integer:: i
      !en = 0.0_wp
      dist = ACO2
      do i = 1,n_CO2
        rl1 = CO2_xyz(i,2,:) - CO2_xyz(i,1,:)
        call pbc(rl1)
        r1sq = sqrt(dot_product(rl1,rl1))
        df1 = KCO2*(r1sq - dist)*rl1(1:NDIM)/r1sq
        CO2_fxyz(i,2,1:NDIM) = CO2_fxyz(i,2,1:NDIM) - df1
        CO2_fxyz(i,1,1:NDIM) = CO2_fxyz(i,1,1:NDIM) + df1
        en(i) = (r1sq - dist)**2
!
        rl2 = CO2_xyz(i,3,:) - CO2_xyz(i,1,:)
        call pbc(rl2)
        r2sq = sqrt(dot_product(rl2,rl2))
        df2 = KCO2*(r2sq - dist)*rl2(1:NDIM)/r2sq
        CO2_fxyz(i,3,1:NDIM) = CO2_fxyz(i,3,1:NDIM) - df2
        CO2_fxyz(i,1,1:NDIM) = CO2_fxyz(i,1,1:NDIM) + df2
        en(i) = en(i) + (r2sq - dist)**2
      end do
      en(1:n_co2) = en(1:n_co2)*0.5_wp*KCO2
    END SUBROUTINE CO2_stretching


    SUBROUTINE CO2_angle_bending(en)
!     Force routine for CO2 bond angle bending
      real(wp),intent(out):: en(:)
      real(wp):: rl1(NDIM),rl2(NDIM),df1(NDIM),df2(NDIM),df3(NDIM)
      real(wp):: r1sq,r2sq,rl1l2,cos_theta,theta,coef,coef1,coef2,coef3
      integer:: i
      en = 0.0_wp
      do i = 1,n_CO2
         rl1 = CO2_xyz(i,2,1:NDIM) - CO2_xyz(i,1,1:NDIM)
         rl2 = CO2_xyz(i,3,1:NDIM) - CO2_xyz(i,1,1:NDIM)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         if (cos_theta >  1.0_wp) cos_theta =  1.0_wp
         if (cos_theta < -1.0_wp) cos_theta = -1.0_wp
         theta = acos(cos_theta)

         if ( (1.0_wp - cos_theta**2) /= 0.0_wp ) then
         coef = K_CO2*(theta - theta_CO2)*(-1.0_wp/sqrt(1.0_wp - cos_theta**2))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
!         CO2_fxyz(i,1,1:NDIM) = CO2_fxyz(i,1,1:NDIM)-coef*(-(rl1+rl2)/coef2 &
!                           + rl1l2*rl1/coef1 &
!                           + rl1l2*rl2/coef3)
!         CO2_fxyz(i,2,1:NDIM) = CO2_fxyz(i,2,1:NDIM)-coef*(rl2/coef2-rl1l2*rl1/coef1)
!         CO2_fxyz(i,3,1:NDIM) = CO2_fxyz(i,3,1:NDIM)-coef*(rl1/coef2-rl1l2*rl2/coef3)
         df1(1:NDIM) =  coef*( (rl1 + rl2)/coef2 - rl1l2*rl1/coef1 - rl1l2*rl2/coef3 )
         df2(1:NDIM) =  coef*( rl1l2*rl1/coef1 - rl2/coef2 )
         df3(1:NDIM) =  coef*( rl1l2*rl2/coef3 - rl1/coef2 )

         CO2_fxyz(i,1,1:NDIM) = CO2_fxyz(i,1,1:NDIM) + df1
         CO2_fxyz(i,2,1:NDIM) = CO2_fxyz(i,2,1:NDIM) + df2
         CO2_fxyz(i,3,1:NDIM) = CO2_fxyz(i,3,1:NDIM) + df3
         end if
         en(i) = (theta - theta_CO2)**2
      end do
      en = en*0.5_wp*K_CO2
   END SUBROUTINE CO2_angle_bending


!  call GAS_stretching(n_GAS,n_atom_GAS,GAS_xyz,ASiO,kSiO)
   SUBROUTINE GAS_stretching(nmol,nat,GAS_xyz,dist,Kgas,en)
!     Force routine for Keating bond stretching
      real(wp),intent(out):: en(:)
      integer,intent(in):: nmol,nat
      real(wp),intent(in):: GAS_xyz(:,:,:),dist,Kgas
      real(wp):: rl1(NDIM),df1(NDIM),r1sq
      integer:: i,j
!real(wp):: u1,u2,dr1(NDIM),dr2(NDIM),Fdr,ro
      en = 0.0_wp
      do i = 1,nmol
!print *
!print *,'mol = ',i
      do j = 2,nat
         rl1 = GAS_xyz(i,j,:) - GAS_xyz(i,1,:)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         df1 = Kgas*(r1sq - dist)*rl1(1:NDIM)/r1sq
         GAS_fxyz(i,j,1:NDIM) = GAS_fxyz(i,j,1:NDIM) - df1
         GAS_fxyz(i,1,1:NDIM) = GAS_fxyz(i,1,1:NDIM) + df1
         en(i) = en(i) + (r1sq - dist)**2
!write(777,*) r1sq
!ro = r1sq
!u1 = 0.5_wp*Kgas*(r1sq-dist)**2
!dr1 = (/  0.1_wp, 0.5_wp, -2.7_wp /)*(-1.0e-10_wp)
!dr2 = (/ -1.0_wp, 0.1_wp,  0.7_wp /)*(-1.0e-10_wp)
!rl1 = GAS_xyz(i,j,:)+dr1 - (GAS_xyz(i,1,:)+dr2)
!call pbc(rl1)
!r1sq = sqrt(dot_product(rl1,rl1))
!Fdr = dot_product(-df1,dr1) + dot_product(df1,dr2)
!u2 = 0.5_wp*Kgas*(r1sq-dist)**2
!print *
!print *,'bond: ',j,' - ',1
!print *,'r1sq-dist = ',r1sq-dist
!print *,'energy = ',0.5_wp*Kgas*(r1sq-dist)**2
!print *
!print *,'   dU = ',u2-u1
!print *,'-F.dr = ',-Fdr
!print *,'|F||dr| = ',sqrt(dot_product(df1,df1))*(r1sq-ro)
!print *
      end do
      end do
      en(1:nmol) = en(1:nmol)*0.5_wp*Kgas
   END SUBROUTINE GAS_stretching


! call GAS_angle_bending(n_GAS,n_atom_GAS,GAS_xyz,ctheta = cos_OSiO,ktheta = KOSiO)
   SUBROUTINE GAS_angle_bending(nmol,nat,GAS_xyz,ctheta,Ktheta,en)
      real(wp),intent(out):: en(:)
      integer,intent(in):: nmol,nat
      real(wp),intent(in):: GAS_xyz(:,:,:),ctheta,Ktheta
      real(wp):: rl1(NDIM),rl2(NDIM),df2(NDIM),df3(NDIM)
      real(wp):: r1sq,r2sq,rl1l2,cos_theta,coef,coef1,coef2,coef3
      integer:: i,j,k
!real(wp):: u1,u2,dr1(NDIM),dr2(NDIM),dr3(NDIM),Fdr,ro
!
      en = 0.0_wp
      do i = 1, nmol
!print *,'mol = ',i
         do j = 2, nat - 1
         do k = j + 1, nat
!print *
!print *,'angle: ',j,' - ',1,' - ',k
            rl1 = GAS_xyz(i,j,1:NDIM) - GAS_xyz(i,1,1:NDIM)
            rl2 = GAS_xyz(i,k,1:NDIM) - GAS_xyz(i,1,1:NDIM)
            call pbc(rl1)
            call pbc(rl2)
            r1sq = sqrt(dot_product(rl1,rl1))
            r2sq = sqrt(dot_product(rl2,rl2))
            rl1l2 = dot_product(rl1,rl2)
            cos_theta = (rl1l2)/(r1sq*r2sq)
            coef = ktheta*(cos_theta - ctheta)
            coef1 = r1sq**3*r2sq
            coef2 = r1sq*r2sq
            coef3 = r1sq*r2sq**3
            df2 = coef*(rl1l2*rl1/coef1 - rl2/coef2)
            df3 = coef*(rl1l2*rl2/coef3 - rl1/coef2)
            en(i) = en(i) + (cos_theta - ctheta)**2
            GAS_fxyz(i,j,1:NDIM) = GAS_fxyz(i,j,1:NDIM) + df2
            GAS_fxyz(i,k,1:NDIM) = GAS_fxyz(i,k,1:NDIM) + df3
            GAS_fxyz(i,1,1:NDIM) = GAS_fxyz(i,1,1:NDIM) - df2 - df3
!write(777,*) cos_theta,cos_theta-ctheta,acos(cos_theta)*180.0/(4.0*atan(1.0))
!print '(a,3f12.6)','rj - r1',rl1
!print '(a,3f12.6)','rk - r1',rl2
!print '(a,3f12.6)','fj',df2
!print '(a,3f12.6)','fk',df3
!print '(a,3f12.6)','f1', -df2 - df3
!u1 = 0.5_wp*ktheta*(cos_theta-ctheta)**2
!dr1 = (/  0.1_wp, 0.5_wp, -2.7_wp /)*(-1.0e-10_wp)
!dr2 = (/ -1.0_wp, 0.1_wp,  0.7_wp /)*(-1.0e-10_wp)
!dr3 = (/ -1.0_wp, 0.3_wp,  0.1_wp /)*(-1.0e-10_wp)
!rl1 = GAS_xyz(i,j,1:NDIM)+dr2 - (GAS_xyz(i,1,1:NDIM)+dr1)
!rl2 = GAS_xyz(i,k,1:NDIM)+dr3 - (GAS_xyz(i,1,1:NDIM)+dr1)
!call pbc(rl1)
!call pbc(rl2)
!r1sq = sqrt(dot_product(rl1,rl1))
!r2sq = sqrt(dot_product(rl2,rl2))
!rl1l2 = dot_product(rl1,rl2)
!cos_theta = (rl1l2)/(r1sq*r2sq)
!coef = ktheta*(cos_theta-ctheta)
!coef1 = r1sq**3*r2sq
!coef2 = r1sq*r2sq
!coef3 = r1sq*r2sq**3
!df2 = coef*(rl1l2*rl1/coef1 - rl2/coef2)
!df3 = coef*(rl1l2*rl2/coef3 - rl1/coef2)
!u2 = 0.5_wp*ktheta*(cos_theta-ctheta)**2
!Fdr = dot_product(-df2-df3,dr1) + dot_product(df2,dr2) + dot_product(df3,dr3)
!print *,'cos_theta-ctheta = ',cos_theta-ctheta
!print *,'energy = ',0.5_wp*ktheta*(cos_theta-ctheta)**2
!print *
!print *,'   dU = ',u2-u1
!print *,'-F.dr = ',-Fdr
!print *
         end do
         end do
      end do
      en = en*0.5_wp*ktheta
   END SUBROUTINE

   SUBROUTINE FORCE_ENERGY_LJ_EL_GAS(atomfirst,atomlast,sysfirst,syslast,nat_gas,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast,nat_gas
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         image_GAS_xyz(M,1,:) = GAS_xyz(M,1,:)
         rr(1:NDIM) = GAS_xyz(M,2,:) - GAS_xyz(M,1,:)
         call pbc(rr)
         image_GAS_xyz(M,2,:) = GAS_xyz(M,1,:) + rr(:)
         rr(1:NDIM) = GAS_xyz(M,3,:) - GAS_xyz(M,1,:)
         call pbc(rr)
         image_GAS_xyz(M,3,:) = GAS_xyz(M,1,:) + rr(:)

         do n = 1,nat_gas
         rr(1:NDIM) = GAS_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = rr(1)**2 + rr(2)**2 + rr(NDIM)**2
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            dele = rr6*(aij(atom(L),GAS_atom(M,n))*rr6 - bij(atom(L),GAS_atom(M,n))) &
                 - vdwcut(atom(L),GAS_atom(M,n)) - vdwsf(atom(L),GAS_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),GAS_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7) &
               + 12.0_wp*aij(atom(L),GAS_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            GAS_fxyz(M,n,1:NDIM) = GAS_fxyz(M,n,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if

         end do

      end do atom_loop
      end do outer_atom_loop
   END SUBROUTINE FORCE_ENERGY_LJ_EL_GAS


   SUBROUTINE FORCE_ENERGY_LJ(sysfirst,syslast,gasfirst,gaslast,nat_gas, &
                             atom_gas,rxyz_gas,fxyz_gas,ULJEL)
      integer,intent(in):: sysfirst,syslast,gasfirst,gaslast,nat_gas
      integer,intent(in):: atom_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      solid_loop: do L = sysfirst,syslast
      mol_loop: do M = gasfirst,gaslast
      atom_loop: do n = 1,nat_gas
         rr(1:NDIM) = rxyz_gas(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         if (r2 < rcut2) then
            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7
            dele = rr6*(aij(atom(L),atom_gas(M,n))*rr6 - bij(atom(L),atom_gas(M,n))) &
                 -   vdwcut(atom(L),atom_gas(M,n)) - vdwsf(atom(L),atom_gas(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),atom_gas(M,n))*(rr7 - rcut7) &
               + 12.0_wp*aij(atom(L),atom_gas(M,n))*(rcut13 - rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            fxyz_gas(M,n,1:NDIM) = fxyz_gas(M,n,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if
      end do atom_loop
      end do mol_loop
      end do solid_loop
   END SUBROUTINE FORCE_ENERGY_LJ


!   PURE SUBROUTINE FORCE_ENERGY_LJ_ATTRACTIVE(nlist,list,atom,rxyz,fxyz,ULJEL)
   PURE SUBROUTINE FORCE_ENERGY_LJ_LIST(nlist,list,atom,rxyz,fxyz,ULJEL)
      USE coordinates_mod, only: pbc, NDIM
      real(wp),parameter:: rcut = 17.0_wp/Angstrom
      real(wp),parameter:: rc0 = 0.5_wp/Angstrom
      real(wp),parameter:: rcut2 = rcut**2
      real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
      real(wp),parameter:: rcut7 = rcut6/rcut
      real(wp),parameter:: rcut13 = rcut7*rcut7*rcut
      integer,intent(in):: list(:),nlist
      integer,intent(in):: atom(:)
      real(wp),intent(in):: rxyz(:,:)
      real(wp),intent(inout):: fxyz(:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13 ,A_ij,B_ij,eij,rstij,vdwcut_ij,vdwsf_ij
      integer:: i,j,L,M
      ULJEL = 0.0_wp
      do L = 1,nlist-1
      i = list(L)
      do M = L+1,nlist
         j = list(M)
         rr(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         if (r2 < rcut2) then
eij = sqrt(epsLJ(atom(i))*epsLJ(atom(j)))
rstij = (sigLJ(atom(i)) + sigLJ(atom(j)))*0.5
eij = eij*10.0_wp
rstij = 1.0/Angstrom
A_ij = 4.0*(rstij**12)*eij
B_ij = 4.0*(rstij**6)*eij
vdwcut_ij = rcut6*(A_ij*rcut6 - B_ij)
vdwsf_ij = rcut6*(-12.0_wp*A_ij*rcut6 + 6.0*B_ij)/rcut

            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7

            dele = rr6*(A_ij*rr6 - B_ij) &
                 -   vdwcut_ij - vdwsf_ij*(r1 - rcut)
            df = (6.0_wp*B_ij*(rr7 - rcut7) &
               + 12.0_wp*A_ij*(rcut13 - rr13))*rr(1:NDIM)/r1

!            dele = rr6*(aij(atom(i),atom(j))*rr6 - bij(atom(i),atom(j))) &
!                 -   vdwcut(atom(i),atom(j)) - vdwsf(atom(i),atom(j))*(r1 - rcut)
!            df = (6.0_wp*bij(atom(i),atom(j))*(rr7 - rcut7) &
!               + 12.0_wp*aij(atom(i),atom(j))*(rcut13 - rr13))*rr(1:NDIM)/r1

!         dele = rr6*(-bij(atom(i),atom(j))) - rcut6*(-bij(atom(i),atom(j))) &
!              - (rcut6*(6.0*-bij(atom(i),atom(j)))/rcut)*(r1 - rcut)
!         df = (6.0_wp*bij(atom(i),atom(j))*(rr7 - rcut7))*rr(1:ndim)/r1
            fxyz(i,1:NDIM) = fxyz(i,1:NDIM) + df
            fxyz(j,1:NDIM) = fxyz(j,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if
      end do
      end do
   END SUBROUTINE FORCE_ENERGY_LJ_LIST


   PURE SUBROUTINE ENERGY_LJ_LIST(ia,nlist,list,atom,rxyz,ULJEL)
      USE coordinates_mod, only: pbc, NDIM
      real(wp),parameter:: rcut = 17.0_wp/Angstrom
      real(wp),parameter:: rc0 = 0.5_wp/Angstrom
      real(wp),parameter:: rcut2 = rcut**2
      real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
      real(wp),parameter:: rcut7 = rcut6/rcut
      real(wp),parameter:: rcut13 = rcut7*rcut7*rcut
      integer,intent(in):: list(:),nlist,ia
      integer,intent(in):: atom(:)
      real(wp),intent(in):: rxyz(:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),r2,r1,rr6,dele,rr7,rr13 ,A_ij,B_ij,eij,rstij,vdwcut_ij,vdwsf_ij
      integer:: j,M
      ULJEL = 0.0_wp
      do M = 1,nlist
         j = list(M)
         if (j == ia) cycle
         rr(1:NDIM) = rxyz(j,1:NDIM) - rxyz(ia,1:NDIM)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         if (r2 < rcut2) then
eij = sqrt(epsLJ(atom(ia))*epsLJ(atom(j)))
rstij = (sigLJ(atom(ia)) + sigLJ(atom(j)))*0.5
eij = eij*10.0_wp
rstij = 1.0/Angstrom
A_ij = 4.0*(rstij**12)*eij
B_ij = 4.0*(rstij**6)*eij
vdwcut_ij = rcut6*(A_ij*rcut6 - B_ij)
vdwsf_ij = rcut6*(-12.0_wp*A_ij*rcut6 + 6.0*B_ij)/rcut

            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7

            dele = rr6*(A_ij*rr6 - B_ij) &
                 -   vdwcut_ij - vdwsf_ij*(r1 - rcut)

!            dele = rr6*(aij(atom(ia),atom(j))*rr6 - bij(atom(ia),atom(j))) &
!                 -   vdwcut(atom(ia),atom(j)) - vdwsf(atom(ia),atom(j))*(r1 - rcut)

!         dele = rr6*(-bij(atom(ia),atom(j))) - rcut6*(-bij(atom(ia),atom(j))) &
!              - (rcut6*(6.0*-bij(atom(ia),atom(j)))/rcut)*(r1 - rcut)
            ULJEL = ULJEL + dele
         end if
      end do
   END SUBROUTINE ENERGY_LJ_LIST


   PURE SUBROUTINE FORCE_ENERGY_LJ2(ifirst,ilast,atom,rxyz,fxyz,ULJEL)
      USE coordinates_mod, only: pbc, NDIM
      real(wp),parameter:: rcut = 17.0_wp/Angstrom
      real(wp),parameter:: rc0 = 0.5_wp/Angstrom
      real(wp),parameter:: rcut2 = rcut**2
      real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
      real(wp),parameter:: rcut7 = rcut6/rcut
      real(wp),parameter:: rcut13 = rcut7*rcut7*rcut
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: atom(:)
      real(wp),intent(in):: rxyz(:,:)
      real(wp),intent(inout):: fxyz(:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13 ,A_ij,B_ij,eij,rstij,vdwcut_ij,vdwsf_ij
      integer:: i,j
      ULJEL = 0.0_wp
      do i = ifirst,ilast-1
      do j = i+1,ilast
         rr(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         if (r2 < rcut2) then
eij = sqrt(epsLJ(atom(i))*epsLJ(atom(j)))
rstij = (sigLJ(atom(i)) + sigLJ(atom(j)))*0.5
eij = eij*10.0_wp
rstij = 1.0/Angstrom
A_ij = 4.0*(rstij**12)*eij
B_ij = 4.0*(rstij**6)*eij
vdwcut_ij = rcut6*(A_ij*rcut6 - B_ij)
vdwsf_ij = rcut6*(-12.0_wp*A_ij*rcut6 + 6.0*B_ij)/rcut

            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7

            dele = rr6*(A_ij*rr6 - B_ij) &
                 -   vdwcut_ij - vdwsf_ij*(r1 - rcut)
            df = (6.0_wp*B_ij*(rr7 - rcut7) &
               + 12.0_wp*A_ij*(rcut13 - rr13))*rr(1:NDIM)/r1

!            dele = rr6*(aij(atom(i),atom(j))*rr6 - bij(atom(i),atom(j))) &
!                 -   vdwcut(atom(i),atom(j)) - vdwsf(atom(i),atom(j))*(r1 - rcut)
!            df = (6.0_wp*bij(atom(i),atom(j))*(rr7 - rcut7) &
!               + 12.0_wp*aij(atom(i),atom(j))*(rcut13 - rr13))*rr(1:NDIM)/r1

!         dele = rr6*(-bij(atom(i),atom(j))) - rcut6*(-bij(atom(i),atom(j))) &
!              - (rcut6*(6.0*-bij(atom(i),atom(j)))/rcut)*(r1 - rcut)
!         df = (6.0_wp*bij(atom(i),atom(j))*(rr7 - rcut7))*rr(1:ndim)/r1
            fxyz(i,1:NDIM) = fxyz(i,1:NDIM) + df
            fxyz(j,1:NDIM) = fxyz(j,1:NDIM) - df
            ULJEL = ULJEL + dele
         end if
      end do
      end do
   END SUBROUTINE FORCE_ENERGY_LJ2


   SUBROUTINE ENERGY_LJ2(ia,ifirst,ilast,atom,rxyz,ULJEL)
      USE coordinates_mod, only: pbc, NDIM
      real(wp),parameter:: rcut = 17.0_wp/Angstrom
      real(wp),parameter:: rc0 = 0.5_wp/Angstrom
      real(wp),parameter:: rcut2 = rcut**2
      real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
      real(wp),parameter:: rcut7 = rcut6/rcut
      real(wp),parameter:: rcut13 = rcut7*rcut7*rcut
      integer,intent(in):: ia,ifirst,ilast
      integer,intent(in):: atom(:)
      real(wp),intent(in):: rxyz(:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(NDIM),r2,r1,rr6,dele,rr7,rr13 ,A_ij,B_ij,eij,rstij,vdwcut_ij,vdwsf_ij
      integer:: j
      ULJEL = 0.0_wp
      do j = ifirst,ilast
         if (j == ia) cycle
         rr(1:NDIM) = rxyz(j,1:NDIM) - rxyz(ia,1:NDIM)
         call pbc(rr)
         r2 = dot_product(rr,rr)
         if (r2 < rcut2) then
eij = sqrt(epsLJ(atom(ia))*epsLJ(atom(j)))
rstij = (sigLJ(atom(ia)) + sigLJ(atom(j)))*0.5
eij = eij*10.0_wp
rstij = 1.0/Angstrom
A_ij = 4.0*(rstij**12)*eij
B_ij = 4.0*(rstij**6)*eij
vdwcut_ij = rcut6*(A_ij*rcut6 - B_ij)
vdwsf_ij = rcut6*(-12.0_wp*A_ij*rcut6 + 6.0*B_ij)/rcut

            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7

            dele = rr6*(A_ij*rr6 - B_ij) &
                 -   vdwcut_ij - vdwsf_ij*(r1 - rcut)

!            dele = rr6*(aij(atom(ia),atom(j))*rr6 - bij(atom(ia),atom(j))) &
!                 -   vdwcut(atom(ia),atom(j)) - vdwsf(atom(ia),atom(j))*(r1 - rcut)
!         dele = rr6*(-bij(atom(ia),atom(j))) - rcut6*(-bij(atom(ia),atom(j))) &
!              - (rcut6*(6.0*-bij(atom(ia),atom(j)))/rcut)*(r1 - rcut)
            ULJEL = ULJEL + dele
         end if
      end do
   END SUBROUTINE ENERGY_LJ2


   SUBROUTINE FORCE_ENERGY_LJ_NLIST(gasfirst,gaslast,nat_gas, &
              atom_gas,rxyz_gas,fxyz_gas,ULJEL)

      integer,intent(in):: gasfirst,gaslast,nat_gas
      integer,intent(in):: atom_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: ULJEL(:)
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer,parameter:: nlm = 5
      integer:: L,M,n,ic,jj,nc,ncell((2*nlm + 1)**3),neigh,nl
      nl = nlayers_ll(rcut - tiny(1.0))
      if (nl > nlm) STOP 'nl > nlm '
      ULJEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
      atom_loop: do n = 1,nat_gas
         ic = CELL( rxyz_gas(M,n,:) )
if (ic < 1) stop 'FORCE_ENERGY_LJ_NLIST: ic < 1'
if (ic > ncelt) stop 'FORCE_ENERGY_LJ_NLIST: ic > ncelt'
         call NEIGCELL(ic,nl,neigh,ncell)
         cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         L = HOC(nc)
         cell_atom_loop: do while (L /= 0)
            rr(1:NDIM) = rxyz_gas(M,n,:) - rxyz(L,:)
            call pbc(rr)
            r2 = dot_product(rr,rr)
            if (r2 < rcut2) then
            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7
            dele = rr6*(aij(atom(L),atom_gas(M,n))*rr6 - bij(atom(L),atom_gas(M,n))) &
                 -   vdwcut(atom(L),atom_gas(M,n)) - vdwsf(atom(L),atom_gas(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),atom_gas(M,n))*(rr7 - rcut7) &
               + 12.0_wp*aij(atom(L),atom_gas(M,n))*(rcut13 - rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            fxyz_gas(M,n,1:NDIM) = fxyz_gas(M,n,1:NDIM) - df
            ULJEL(M) = ULJEL(M) + dele
            end if
         L = LL(L)
         end do cell_atom_loop
         end do cell_loop
      end do atom_loop
      end do mol_loop
   END SUBROUTINE FORCE_ENERGY_LJ_NLIST

END MODULE Lj_el_mod


!!>include 'force_en_el_mod.f90'

MODULE force_en_el_mod
    USE precision_mod
    USE coordinates_mod
    USE charges_mod
    USE constants_mod
    implicit none

CONTAINS

   SUBROUTINE FORCE_ENERGY_EL(sysfirst,syslast,gasfirst,gaslast,nat_gas, &
                             charge_gas,rxyz_gas,fxyz_gas,UEL)
      integer,intent(in):: sysfirst,syslast,gasfirst,gaslast,nat_gas
      real(wp),intent(in):: charge_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: UEL
      real(wp):: rgas(NDIM),rr(NDIM),df(NDIM),r2,r1,du
      integer:: L,M,n
      UEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
         atom_loop: do n = 1,nat_gas
            rr = rxyz_gas(M,n,:) - rxyz_gas(M,1,:)
            call pbc(rr)
            rgas = rxyz_gas(M,1,:) + rr
            solid_loop: do L = sysfirst,syslast
               rr(1:NDIM) = rgas - rxyz(L,:)
               call pbc(rr)
               r2 = dot_product(rr,rr)
               r1 = sqrt(r2)
               du = charge(L)*charge_gas(M,n)/r1
               df = du*rr/r2
               fxyz(L,1:NDIM) = fxyz(L,1:NDIM) - df
               fxyz_gas(M,n,1:NDIM) = fxyz_gas(M,n,1:NDIM) + df
               Uel = Uel + du
            end do solid_loop
         end do atom_loop
      end do mol_loop
   END SUBROUTINE


   SUBROUTINE FORCE_ENERGY_EL_X2_dummy(sysfirst,syslast,gasfirst,gaslast, &
                                       charge_gas,rxyz_gas,fxyz_gas,UEL)
      integer,intent(in):: sysfirst,syslast,gasfirst,gaslast
      real(wp),intent(in):: charge_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: UEL
      real(wp):: rgas(3,NDIM),rr(NDIM),df(NDIM),r2,r1,du
      integer:: L,M,n
      UEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
         rgas(1,:) = rxyz_gas(M,1,:)
         rr(1:NDIM) = rxyz_gas(M,2,:) - rxyz_gas(M,1,:)
         call pbc(rr)
         rgas(2,:) = rxyz_gas(M,1,:) + rr(:)
         rgas(3,:) = (rgas(1,:) + rgas(2,:))*0.5_wp
         solid_loop: do L = sysfirst,syslast
            do n = 1,2
               rr(1:NDIM) = rgas(n,:) - rxyz(L,:)
               call pbc(rr)
               r2 = dot_product(rr,rr)
               r1 = sqrt(r2)
               du = charge(L)*charge_gas(M,n)/r1
               df = du*rr/r2
               fxyz(L,1:NDIM) = fxyz(L,1:NDIM) - df
               fxyz_gas(M,n,1:NDIM) = fxyz_gas(M,n,1:NDIM) + df
               Uel = Uel + du
            end do
            rr(1:NDIM) = rgas(3,:) - rxyz(L,:)
            call pbc(rr)
            r2 = dot_product(rr,rr)
            r1 = sqrt(r2)
            du = charge(L)*charge_gas(M,3)/r1
            df = du*rr/r2
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) - df
            fxyz_gas(M,1,1:NDIM) = fxyz_gas(M,1,1:NDIM) + 0.5_wp*df
            fxyz_gas(M,2,1:NDIM) = fxyz_gas(M,2,1:NDIM) + 0.5_wp*df
            UEL = Uel + du
         end do solid_loop
      end do mol_loop
   END SUBROUTINE FORCE_ENERGY_EL_X2_dummy

END MODULE

!!>include 'tcf_mod.f90'

MODULE tcf_mod
   USE precision_mod, only: wp
   implicit none
!
   TYPE tcf_type
      integer:: icount
      integer:: isamp
      integer:: itau
      integer:: nfreq
      integer:: ntau
      integer:: nmol
      integer:: iout
      real(wp),allocatable:: tcfmol(:,:)
      real(wp),allocatable:: hist(:,:,:)
   END TYPE
!
   public:: new_tcf
!
CONTAINS

   SUBROUTINE new_tcf(nmol,ntau,nfreq,T,fname)
      USE files_mod
      TYPE(tcf_type),intent(out):: T
      integer,intent(in):: nmol,ntau,nfreq
      character(*),intent(in):: fname
      T%ntau = ntau
      T%nmol = nmol
      allocate( T%hist(nmol,ntau,3), T%tcfmol(nmol,ntau) )
      T%hist = 0.0
      T%tcfmol = 0.0
      T%icount = 0
      T%isamp = 0
      T%itau = 0
      T%nfreq = nfreq
      call get_free_file_unit(T%iout)
      open(unit=T%iout, file = trim(adjustl(fname)))
   END SUBROUTINE new_tcf

END MODULE tcf_mod


!!>include 'msd_mod3.f90'

MODULE msd_mod
   USE precision_mod, only: wp
   USE tcf_mod
   implicit none
   public:: msd_samp, msd_print
!
CONTAINS

   SUBROUTINE msd_samp(T,xyz)
      real(wp),intent(in):: xyz(:,:)
      type(tcf_type),intent(inout):: T
      integer:: i,L,m,nt
!
      T%icount = T%icount + 1
      T%itau = T%itau + 1
      if (T%itau > T%ntau) T%itau = T%itau - T%ntau
      nt = T%itau
! print*,'T%icount =',T%icount 
! print*,'T%itau =',T%itau 
 print*,'nt =',nt 
!
      do i = 1,T%nmol
         T%hist(i,nt,:) = xyz(i,:)
      end do
!
!print*,',T%nmol',T%nmol
!print*,'T%ntau =',T%ntau
      if (T%icount <= T%ntau) RETURN
      T%isamp = T%isamp + 1
print*,'T%isamp======',T%isamp 
!if (T%isamp == 5)stop 'T%isamp == 5'
      do L = 1,T%ntau
!print*,'L==',L      
         m = nt - L + 1
print*,'m==',m
         if (m <= 0) m = m + T%ntau
print*,'m=======',m         
         do i = 1,T%nmol
            T%tcfmol(i,L) = T%tcfmol(i,L) &
                          + ( T%hist(i,nt,1) - T%hist(i,m,1) )**2 &
                          + ( T%hist(i,nt,2) - T%hist(i,m,2) )**2 &
                          + ( T%hist(i,nt,3) - T%hist(i,m,3) )**2
!print*,'T%tcfmol(i,L) =', T%tcfmol(i,L)                          
         end do
      end do
   END SUBROUTINE msd_samp


   SUBROUTINE msd_print(T, dtfs)
      type(tcf_type),intent(in):: T
      real(wp),intent(in):: dtfs
      real(wp):: msd(size(T%tcfmol,dim=2))
      real(wp):: rnmal,tau
      integer:: j
      if (T%isamp <= 0) RETURN
      rnmal = 1.0_wp/real(T%isamp,wp)
! sum individual molecule tracer msds
      msd(:) = 0.0_wp
      do j = 1,T%nmol
         msd(:) = msd(:) + T%tcfmol(j,:)
      end do
      msd(:) = msd(:)*rnmal/T%nmol
!
      do j = 1, T%ntau
         tau = (j - 1)*T%nfreq*dtfs*0.001_wp
         write(T%iout,*) tau,msd(j)
      end do
      write(T%iout,'(/)')
   END SUBROUTINE msd_print

   SUBROUTINE msd_sample(T,xyz)
      real(wp),intent(in):: xyz(:,:)
      type(tcf_type),intent(inout):: T
      call msd_samp_arr(T%nmol, T%ntau, T%itau, T%icount, T%isamp, xyz, T%tcfmol, T%hist)
   END SUBROUTINE msd_sample

   SUBROUTINE msd_samp_arr(nmol, ntau, itau, icount, isamp, xyz, tcfmol, hist)
      integer,intent(in):: nmol, ntau
      integer,intent(inout):: itau, icount, isamp
      real(wp),intent(in):: xyz(:,:)
      real(wp),intent(inout):: tcfmol(:,:),hist(:,:,:)
      integer:: i,L,m
!
      icount = icount + 1
      itau = itau + 1
      if (itau > ntau) itau = itau - ntau
!
      do i = 1,nmol
         hist(i,itau,:) = xyz(i,:)
      end do
!
      if (icount <= ntau) RETURN
      isamp = isamp + 1
      do L = 1,ntau
         m = itau - L + 1
         if (m <= 0) m = m + ntau
         do i = 1,nmol
            tcfmol(i,L) = tcfmol(i,L) &
                        + ( hist(i,itau,1) - hist(i,m,1) )**2 &
                        + ( hist(i,itau,2) - hist(i,m,2) )**2 &
                        + ( hist(i,itau,3) - hist(i,m,3) )**2
         end do
      end do
   END SUBROUTINE msd_samp_arr

END MODULE msd_mod


!!>include 'vcf_mod3.f90'

MODULE vcf_mod
   USE precision_mod, only: wp
   USE tcf_mod
   implicit none
   public:: vcf_samp, vcf_print
!
CONTAINS

   SUBROUTINE vcf_samp(T,vcm)
      real(wp),intent(in):: vcm(:,:)
      type(tcf_type),intent(inout):: T
      integer:: i,L,m,nt
!
      T%icount = T%icount + 1
      T%itau = T%itau + 1
      if (T%itau > T%ntau) T%itau = T%itau - T%ntau
      nt = T%itau
!
      do i = 1,T%nmol
         T%hist(i,nt,:) = vcm(i,:)
      end do
!
      if (T%icount <= T%ntau) RETURN
      T%isamp = T%isamp + 1
      do L = 1,T%ntau
         m = nt - L + 1
         if (m <= 0) m = m + T%ntau
         do i = 1,T%nmol
            T%tcfmol(i,L) = T%tcfmol(i,L) &
                          + T%hist(i,nt,1)*T%hist(i,m,1) &
                          + T%hist(i,nt,2)*T%hist(i,m,2) &
                          + T%hist(i,nt,3)*T%hist(i,m,3)
         end do
      end do
   END SUBROUTINE vcf_samp


   SUBROUTINE vcf_print(T, dtfs)
      type(tcf_type),intent(in):: T
      real(wp),intent(in):: dtfs
      real(wp):: vcf(size(T%tcfmol,dim=2))
      real(wp):: rnmal,tau
      integer:: j
      if (T%isamp <= 0) RETURN
      rnmal = 1.0_wp/real(T%isamp,wp)
! sum individual molecule tracer VCFs
      vcf(:) = 0.0_wp
      do j = 1,T%nmol
         vcf(:) = vcf(:) + T%tcfmol(j,:)
      end do
!
      write(T%iout,"('# ',g15.7)") vcf(1)*rnmal/T%nmol
      do j = 1, T%ntau
         tau = (j - 1)*T%nfreq*dtfs*0.001_wp
         write(T%iout,*) tau,vcf(j)/vcf(1)
      end do
      write(T%iout,'(/)')
   END SUBROUTINE vcf_print


   SUBROUTINE vcf_sample(T,vcm)
      real(wp),intent(in):: vcm(:,:)
      type(tcf_type),intent(inout):: T
      call vcf_samp_arr(T%nmol, T%ntau, T%itau, T%icount, T%isamp, vcm, T%tcfmol, T%hist)
   END SUBROUTINE vcf_sample


   SUBROUTINE vcf_samp_arr(nmol, ntau, itau, icount, isamp, vcm, tcfmol, hist)
      integer,intent(in):: nmol, ntau
      integer,intent(inout):: itau, icount, isamp
      real(wp),intent(in):: vcm(:,:)
      real(wp),intent(inout):: tcfmol(:,:),hist(:,:,:)
      integer:: i,L,m
!
      icount = icount + 1
      itau = itau + 1
      if (itau > ntau) itau = itau - ntau
!
      do i = 1,nmol
         hist(i,itau,:) = vcm(i,:)
      end do
!
      if (icount <= ntau) RETURN
      isamp = isamp + 1
      do L = 1,ntau
         m = itau - L + 1
         if (m <= 0) m = m + ntau
         do i = 1,nmol
            tcfmol(i,L) = tcfmol(i,L) &
                        + hist(i,itau,1)*hist(i,m,1) &
                        + hist(i,itau,2)*hist(i,m,2) &
                        + hist(i,itau,3)*hist(i,m,3)
         end do
      end do
   END SUBROUTINE vcf_samp_arr

END MODULE vcf_mod


!!>include 'rdf_mod.f90'

MODULE RDF_MOD
   USE precision_mod, only: wp
   USE global_vars_mod
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nbin = 50
!  real(wp),private:: zl0 = 2.0_wp/Angstrom,zu0 = 22.0_wp/Angstrom
   real(wp),private::  zl = 10.0_wp/Angstrom, zu = 15.0_wp/Angstrom
   real(wp),parameter,private:: rc = 7.0_wp/Angstrom, rc2 = rc**2
   real(wp),private:: RDF(nbin),RDF1(nbin)
   real(wp),private:: RDFT(nbin,3)
   integer,private:: npairs(3),na
   real(wp),private:: DL1,rho,npairstot
   integer,private:: ngr = 1
   PUBLIC:: RDF_CALC, RDF_PRINT
CONTAINS

   SUBROUTINE RDF_CALC0(ifirst,ilast)
      integer,intent(in):: ifirst,ilast
      real(wp):: rv(NDIM),r,r2,rmax,rmax2
      integer:: i,j,k
      rmax = maxval(boxl)*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDF = 0.0_wp
      do i = ifirst,ilast - 1
         do j = i + 1,ilast
            rv(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rmax2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1) + 1
            RDF(k) = RDF(k) + 2.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC0

   SUBROUTINE RDF_CALC(ifirst,ilast)
      integer,intent(in):: ifirst,ilast
      real(wp):: rv(NDIM),r,r2,zi,zj
      integer:: i,j,k
      DL1 = rc/real(nbin,wp)
      RDF = 0.0_wp
      npairstot = 0.0_wp
      na = 0.0_wp
      do i = ifirst,ilast
         zi = rxyz(i,3)
         if ( (zi > zu) .or. (zi < zl) ) cycle
         na = na + 1
         do j = ifirst,ilast
            if (i == j) cycle
            zj = rxyz(j,3)
            if ( (zj > zu + rc) .or. (zj < zl - rc) ) cycle
!            if ( (zi > zu) .or. (zi < zl) ) then
!            if ( (zj > zu) .or. (zj < zl) ) then
!               cycle
!            end if
!            end if
            rv(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rv)
            r2 = dot_product(rv,rv)
!print *,i,j,'###################### ',sqrt(r2),rc
            if (r2 > rc2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write(*,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
               stop
            end if
            k = int(r/DL1) + 1
            RDF(k) = RDF(k) + 1.0_wp
            npairstot = npairstot + 1.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC

   SUBROUTINE RDF_CALC_ATOMTYPES(ifirst,ilast)
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      real(wp):: rv(NDIM),r,r2,rmax,rmax2,zi,zj
      integer:: i,j,k,it,ia,ja
      rmax = maxval(boxl)*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rc/real(nbin,wp)
      RDFT = 0.0_wp
      npairs = 0
      do i = ifirst,ilast
         ia = atom(i)
         zi = rxyz(i,3)
         do j = i + 1,ilast
            ja = atom(j)
            zj = rxyz(j,3)
            if ( (zi > zu) .or. (zi < zl) ) then
            if ( (zj > zu) .or. (zj < zl) ) then
               cycle
            end if
            end if
            if ((ia == iSilicon).and.(ja == iSilicon)) then
               it = 1
            else if ((ia == iOxygen).and.(ja == iOxygen)) then
               it = 2
            else if ((ia == iOxygenH).and.(ja == iOxygen)) then
               it = 2
            else if ((ia == iOxygen).and.(ja == iOxygenH)) then
               it = 2
            else if ((ia == iOxygenH).and.(ja == iOxygenH)) then
               it = 2
            else if ((ia == iSilicon).and.(ja == iOxygen)) then
               it = 3
            else if ((ia == iSilicon).and.(ja == iOxygenH)) then
               it = 3
            else if ((ia == iOxygen).and.(ja == iSilicon)) then
               it = 3
            else if ((ia == iOxygenH).and.(ja == iSilicon)) then
               it = 3
            end if
            rv(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rc2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1) + 1
            RDFT(k,it) = RDFT(k,it) + 1.0_wp
            npairs(it) = npairs(it) + 1
         end do
      end do
   END SUBROUTINE

   SUBROUTINE RDF_CALC_ATOMTYPES0(ifirst,ilast)
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      real(wp):: rv(NDIM),r,r2,rmax,rmax2
      integer:: i,j,k,it,ia,ja
      rmax = maxval(boxl)*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDFT = 0.0_wp
      npairs = 0
      do i = ifirst,ilast - 1
         ia = atom(i)
         do j = i + 1,ilast
            ja = atom(j)
            if ((ia == iSilicon).and.(ja == iSilicon)) then
               it = 1
            else if ((ia == iOxygen).and.(ja == iOxygen)) then
               it = 2
            else if ((ia == iOxygenH).and.(ja == iOxygen)) then
               it = 2
            else if ((ia == iOxygen).and.(ja == iOxygenH)) then
               it = 2
            else if ((ia == iOxygenH).and.(ja == iOxygenH)) then
               it = 2
            else if ((ia == iSilicon).and.(ja == iOxygen)) then
               it = 3
            else if ((ia == iSilicon).and.(ja == iOxygenH)) then
               it = 3
            else if ((ia == iOxygen).and.(ja == iSilicon)) then
               it = 3
            else if ((ia == iOxygenH).and.(ja == iSilicon)) then
               it = 3
            end if
            rv(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rmax2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1) + 1
            RDFT(k,it) = RDFT(k,it) + 1.0_wp
            npairs(it) = npairs(it) + 1
         end do
      end do
   END SUBROUTINE

   SUBROUTINE RDF_PRINT_ATOMTYPES0(io,nattached,ifirst,ilast)
      integer,intent(in):: io,nattached
      integer,intent(in):: ifirst,ilast
      integer:: i,ityp
      write(io,*)'# nattached = ',nattached
      write(io,*)'# natoms = ',ilast-ifirst+1
      write(io,*)'# Si - Si'
      ityp = 1
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i - 0.5_wp)*DL1,RDFT(i,ityp)/npairs(ityp)
      end do
      write(io,'(/)')
      write(io,*)'# O - O'
      ityp = 2
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i - 0.5_wp)*DL1,RDFT(i,ityp)/npairs(ityp)
      end do
      write(io,'(/)')
      write(io,*)'# Si - O'
      ityp = 3
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i - 0.5_wp)*DL1,RDFT(i,ityp)/npairs(ityp)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE

   SUBROUTINE RDF_PRINT_ATOMTYPES(io,nattached,ifirst,ilast)
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: io,nattached
      integer:: i,ityp,nati(2),nat
      real(wp):: r,dv,zi,vol
      nat = 0
      nati = 0
      do i = ifirst,ilast
         zi = rxyz(i,3)
         if ( (zi <= zu) .and. (zi >= zl) ) then
            nat = nat + 1
            if (atom(i) == iSilicon) then
               nati(1) = nati(1) + 1
            else if (atom(i) == iOxygen) then
               nati(2) = nati(2) + 1
            else if (atom(i) == iOxygenH) then
               nati(2) = nati(2) + 1
            end if
         end if
      end do
      vol = ((zu - zl)*boxl(1)*boxl(2))
      write(io,*)'# nattached = ',nattached
      write(io,*)'# natoms = ',ilast-ifirst+1
      write(io,*)'# Si - Si'
      ityp = 1
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdft(i,ityp) = rdft(i,ityp)/(ngr*dv*(npairs(ityp)/vol))
         write(14,'(f12.6,2x,g16.8)') r,RDFT(i,ityp)
      end do
      write(io,'(/)')
      write(io,*)'# O - O'
      ityp = 2
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdft(i,ityp) = rdft(i,ityp)/(ngr*dv*(npairs(ityp)/vol))
         write(14,'(f12.6,2x,g16.8)') r,RDFT(i,ityp)
      end do
      write(io,'(/)')
      write(io,*)'# Si - O'
      ityp = 3
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdft(i,ityp) = rdft(i,ityp)/(ngr*dv*(npairs(ityp)/vol))
         write(14,'(f12.6,2x,g16.8)') r,RDFT(i,ityp)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE

   SUBROUTINE RDF_PRINT(io,nattached,ifirst,ilast)
      integer,intent(in):: io,nattached
      integer,intent(in):: ifirst,ilast
      integer:: i,nat
      real(wp):: r,dv,zi,vol
      nat = 0
      do i = ifirst,ilast
         zi = rxyz(i,3)
         if ( (zi <= zu + rc) .and. (zi >= zl - rc) ) then
            nat = nat + 1
         end if
      end do
print * , 'nat = ',nat
print * , 'na  = ',na
      write(io,*)'# nattached = ',nattached
      write(io,*)'# natoms = ',ilast-ifirst+1
      write(io,*)'# nat = ',nat
      write(io,*)'# npairstot = ',npairstot
      write(io,*)'# npairstot/(nat - 1) = ',npairstot/(nat - 1)
      rho = nat/((zu - zl)*boxl(1)*boxl(2))
      vol = ((zu - zl)*boxl(1)*boxl(2))
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)/nat
      end do
      write(io,'(/)')
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdf(i) = rdf(i)/(ngr*dv*rho*na)
!         rdf(i) = rdf(i)/(ngr*dv*(nat/vol))
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT

   SUBROUTINE RDF_PRINT0(io,nattached,ifirst,ilast)
      integer,intent(in):: ifirst,ilast
      integer,intent(in):: io,nattached
      integer:: i,nat
      real(wp):: r,dv
      nat = ilast-ifirst+1
      write(io,*)'# nattached = ',nattached
      write(io,*)'# natoms = ',nat
      rho = nat/(maxval(rxyz(ifirst:ilast,3))*boxl(1)*boxl(2))
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)/nat
      end do
      write(io,'(/)')
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdf(i) = rdf(i)/(ngr*dv*rho*nat)
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT0

END MODULE RDF_MOD

!!>include 'utility_mod.f90'

MODULE utility_mod
   USE precision_mod
   implicit none

   interface kinetic_energy
      module procedure kinetic_energy_arr2, kinetic_energy_arr3
   end interface kinetic_energy

   interface rem_com_momentum
      module procedure rem_com_momentum_arr2, rem_com_momentum_arr3
   end interface rem_com_momentum

CONTAINS

   PURE FUNCTION kinetic_energy_arr2(nat,vxyz,mass)
         real(wp):: kinetic_energy_arr2
         integer,intent(in):: nat
         real(wp),intent(in):: vxyz(:,:),mass(:,:)
         integer:: j
         kinetic_energy_arr2 = 0.0_wp
         do j = 1,nat
            kinetic_energy_arr2 = kinetic_energy_arr2 + mass(j,1)*dot_product(vxyz(j,:),vxyz(j,:))
         end do
         kinetic_energy_arr2 = kinetic_energy_arr2*0.5_wp
   END FUNCTION kinetic_energy_arr2

   PURE FUNCTION kinetic_energy_arr3(nmol,nat,vxyz,mass)
         real(wp):: kinetic_energy_arr3
         integer,intent(in):: nmol,nat
         real(wp),intent(in):: vxyz(:,:,:),mass(:,:,:)
         integer:: i,j
         kinetic_energy_arr3 = 0.0_wp
         do j = 1,nat
         do i = 1, nmol
            kinetic_energy_arr3 = kinetic_energy_arr3 + mass(i,j,1)*dot_product(vxyz(i,j,:),vxyz(i,j,:))
!           kinetic_energy = kinetic_energy + mass(i,j,1) &
!                          *(vxyz(i,j,1)*vxyz(i,j,1) + vxyz(i,j,2)*vxyz(i,j,2) + vxyz(i,j,3)*vxyz(i,j,3))
         end do
         end do
         kinetic_energy_arr3 = kinetic_energy_arr3*0.5_wp
   END FUNCTION kinetic_energy_arr3


   SUBROUTINE rem_com_momentum_arr2(natom,vxyz,mass)
      integer,intent(in):: natom
      real(wp),intent(in):: mass(:,:)
      real(wp),intent(inout):: vxyz(:,:)
      real(wp):: SUMX,SUMY,SUMZ,sum_mass
      integer:: i
      if (natom <= 0) RETURN
! REMOVE NET MOMENTUM
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      sum_mass = 0.0_wp
      do i = 1, natom
         SUMX = SUMX + mass(i,1)*vxyz(i,1)
         SUMY = SUMY + mass(i,2)*vxyz(i,2)
         SUMZ = SUMZ + mass(i,3)*vxyz(i,3)
         sum_mass = sum_mass + mass(i,1)
      end do
      SUMX = SUMX/sum_mass
      SUMY = SUMY/sum_mass
      SUMZ = SUMZ/sum_mass
      do i = 1, natom
         vxyz(i,1) = vxyz(i,1) - SUMX
         vxyz(i,2) = vxyz(i,2) - SUMY
         vxyz(i,3) = vxyz(i,3) - SUMZ
      end do
   END SUBROUTINE rem_com_momentum_arr2


   SUBROUTINE rem_com_momentum_arr3(nmol,nat,vxyz,mass)
      integer,intent(in):: nmol,nat
      real(wp),intent(in):: mass(:,:,:)
      real(wp),intent(inout):: vxyz(:,:,:)
      real(wp):: SUMX,SUMY,SUMZ,sum_mass
      integer:: i,j
      if (nmol <= 0 .or. nat <= 0) RETURN
! REMOVE NET MOMENTUM
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      sum_mass = 0.0_wp
      do i = 1,nmol
      do j = 1,nat
         SUMX = SUMX + vxyz(i,j,1)*mass(i,j,1)
         SUMY = SUMY + vxyz(i,j,2)*mass(i,j,2)
         SUMZ = SUMZ + vxyz(i,j,3)*mass(i,j,3)
         sum_mass = sum_mass + mass(i,j,1)
      end do
      end do
      SUMX = SUMX/sum_mass
      SUMY = SUMY/sum_mass
      SUMZ = SUMZ/sum_mass
      do i = 1,nmol
      do j = 1,nat
         vxyz(i,j,1) = vxyz(i,j,1) - SUMX
         vxyz(i,j,2) = vxyz(i,j,2) - SUMY
         vxyz(i,j,3) = vxyz(i,j,3) - SUMZ
      end do
      end do
   END SUBROUTINE rem_com_momentum_arr3


!   SUBROUTINE tetra_coords5(bondl,tetra)
!      real(wp),intent(in):: bondl
!      real(wp),intent(out):: tetra(3,5)
!!
!! returns a set of coordinates for a tetrahedral molecule
!! with the central atom at the origin (atom 1), one atom
!! on the z axis below it (atom 2), and atoms 3,4, and 5 on
!! a plane above the origin.
!!
!      real(wp):: zz,sth,cth,rxy,xx,yy
!      zz = (1.0_wp/3.0_wp)*bondl
!      sth = sqrt(3.0_wp)*0.5_wp
!      cth = -0.5_wp
!      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
!      xx = 0.0_wp
!      yy = -rxy
!      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
!      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
!      tetra(1:3,3) = (/ xx, yy, zz /)
!      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
!      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
!   END SUBROUTINE tetra_coords5

END MODULE utility_mod


!!>include 'energy_keating_mod_sio2only.f90'

MODULE energy_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none

CONTAINS

   SUBROUTINE energy_keating_bond(nbond,ibond,kbond,abond,en)
! Calculates the Keating Bond Energy
      USE coordinates_mod
      integer,intent(in):: nbond,ibond(:,:)
      real(wp),intent(in):: kbond(:),abond(:)
      real(wp),intent(out):: en
      integer:: ib
      real(wp):: r1(NDIM),r1sq
!
      en = 0.0_wp
      do ib = 1,nbond  ! Bond Stretch
         r1 = rxyz(ibond(1,ib),1:NDIM) - rxyz(ibond(2,ib),1:NDIM)
         call pbc(r1)
         r1sq = sqrt(dot_product(r1,r1))
         en = en + 0.5_wp*kbond(ib)*(r1sq - Abond(ib))**2
      end do
   END SUBROUTINE energy_keating_bond

   SUBROUTINE energy_keating_angle(nang,iang,kang,cang,en)
! Calculates the Keating Angle Energy
      USE coordinates_mod
      integer,intent(in):: nang,iang(:,:)
      real(wp),intent(in):: kang(:),cang(:)
      real(wp),intent(out):: en
      integer:: ia
      real(wp):: r1(NDIM),r2(NDIM),cos_theta
!
      en = 0.0_wp
      do ia = 1,nang  ! Bond Angle
         r1 = rxyz(iang(1,ia),1:NDIM) - rxyz(iang(2,ia),1:NDIM)
         r2 = rxyz(iang(3,ia),1:NDIM) - rxyz(iang(2,ia),1:NDIM)
         call pbc(r1)
         call pbc(r2)
         cos_theta = dot_product(r1,r2)/(sqrt(dot_product(r1,r1))*sqrt(dot_product(r2,r2)))
         en = en + 0.5_wp*kang(ia)*(cos_theta - cang(ia))**2
      end do
   END SUBROUTINE energy_keating_angle

   SUBROUTINE energy_keating(en)
! Calculates the total Keating Energy
! by summing over the bond & angle lists
      USE bond_angle_list_mod
      USE coordinates_mod
      real(wp),intent(out):: en
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: r1(NDIM),r1sq,r2(NDIM)
      real(wp):: cos_theta
!
      en = 0.0_wp
      do ib = 1,nbondtot  ! Bond Stretch
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         r1sq = sqrt(dot_product(r1,r1))
         en = en + 0.5_wp*kbond(ib)*(r1sq-Abond(ib))**2
      end do
!
      do ia = 1,nang  ! Bond Angle
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         r2 = rxyz(ia3,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         call pbc(r2)
         cos_theta = dot_product(r1,r2)/ &
            (sqrt(dot_product(r1,r1))*sqrt(dot_product(r2,r2)))
         en = en + 0.5_wp*kang(ia)*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE energy_keating


   PURE FUNCTION energy_bond(r1,r2,d0,kb)
      USE coordinates_mod
      real(wp):: energy_bond
      real(wp),intent(in):: r1(:),r2(:),d0,kb
      real(wp):: r12(NDIM),r12l
      r12 = r1 - r2
      call pbc(r12)
!     r12l = sqrt(dot_product(r12,r12))
      r12l = sqrt( r12(1)**2 + r12(2)**2 + r12(3)**2 )
      energy_bond = 0.5_wp*kb*(r12l - d0)**2
   END FUNCTION energy_bond


   PURE FUNCTION energy_angle(r1,r2,r3,c0,ka)
      USE coordinates_mod
      real(wp):: energy_angle
      real(wp),intent(in):: r1(:),r2(:),r3(:),c0,ka
      real(wp):: r21(NDIM),r23(NDIM),cos_theta
      r21 = r1 - r2
      r23 = r3 - r2
      call pbc(r21)
      call pbc(r23)
!     cos_theta = dot_product(r21,r23)/ &
!          (sqrt(dot_product(r21,r21))*sqrt(dot_product(r23,r23)))
      cos_theta = (r21(1)*r23(1)+r21(2)*r23(2)+r21(3)*r23(3))/ &
         (sqrt(r21(1)*r21(1)+r21(2)*r21(2)+r21(3)*r21(3))* &
          sqrt(r23(1)*r23(1)+r23(2)*r23(2)+r23(3)*r23(3)))
      energy_angle = 0.5_wp*ka*(cos_theta - c0)**2
   END FUNCTION energy_angle


!  PURE
   FUNCTION keating_energy_atom(ia) result(en)
! Calculates the Keating energy of all the bonds & angles containing atom ia
      USE atom_types_mod
      USE coordinates_mod
      USE connectivity_mod
      USE bond_angle_types_mod
!     USE bond_angle_types_mod
      real(wp):: en
      integer,intent(in):: ia
      integer:: j,jj,k,kk
      real(wp):: c0,ck !,kb,d0
!
      en = 0.0_wp
      select case(atom(ia))
      case(iSilicon)
         ck = KOSiO
         c0 = cos_OSiO
      case(iOxygen)
         ck = KSiOSi
         c0 = cos_SiOSi
      case default
      end select
! the angles with ia as centre atom (j--ia--k)
      do jj = 1,ncmax(atom(ia))-1
         j = proximity(ia,jj)
         if (j == 0) cycle
         do kk = jj + 1,ncmax(atom(ia))
            k = proximity(ia,kk)
            if (k == 0) cycle
            ! call angle_type(atom(j),atom(ia),atom(k),c0,ck)
            en = en + energy_angle( rxyz(j,:),rxyz(ia,:),rxyz(k,:),c0,ck )
         end do
      end do

! bonds and other angles involving atom ia
      do jj = 1,ncmax(atom(ia))
         j = proximity(ia,jj)
         if ( j == 0) CYCLE
         ! call bond_type(atom(ia),atom(j),d0,kb)
         en = en + energy_bond( rxyz(ia,:),rxyz(j,:),ASiO,KSiO)
         select case(atom(j))
         case(iSilicon)
            ck = KOSiO
            c0 = cos_OSiO
         case(iOxygen)
            ck = KSiOSi
            c0 = cos_SiOSi
         case default
         end select
         do kk = 1,ncmax(atom(j))
            k = proximity(j,kk)
            if (k == 0) CYCLE
            if (k == ia) CYCLE
            ! call angle_type(atom(ia),atom(j),atom(k),c0,ck)
            en = en + energy_angle( rxyz(ia,:),rxyz(j,:),rxyz(k,:),c0,ck )
         end do
      end do
   END FUNCTION keating_energy_atom


!  PURE
   FUNCTION keating_energy_List(nl,list) result(en)
! Calculate the total bond & angle energy involving atoms 'i' in list
! interacting with atoms j & k either inside or outside the list
      USE atom_types_mod
      USE coordinates_mod
      USE connectivity_mod
      USE bond_angle_types_mod
      USE Keating_parameters_mod
      USE list_mod
      real(wp):: en
      integer,intent(in):: nl,list(:)
      real(wp):: ck,c0 !,d0,kb
      real(wp):: eang,ebnd
      integer:: i,il,j,jj,k,kk
!
      en = 0.0_wp
      list_loop: do il = 1, nl
      i = list(il)

      ! loop over atoms j bonded to i
      do jj = 1, ncmax(atom(i))

         j = proximity(i,jj)
         if (j == 0) cycle
         ! bond i-j
         ! call bond_type(atom(i),atom(j),d0,kb)
         ebnd = energy_bond(rxyz(i,:),rxyz(j,:),ASiO,KSiO)
         en = en + 0.5_wp*ebnd

         ! if (ANY(list(1:nl) == j)) cycle
         if (in_list(j,list,nl)) cycle
         !
         ! if j is not in list
         !          |     k
         !  list    |   /
         !       i--|--j----k
         !          |   \
         !          |    k  outside list
         !
         ! then must include angles i-j-k
         ! and add the remaining 1/2 of the bond energy
         !
         en = en + 0.5_wp*ebnd
         ! call angle_type(atom(i),atom(j),atom(k),c0,ck)
         select case(atom(j))
         case(iSilicon)
            ck = KOSiO
            c0 = cos_OSiO
         case(iOxygen)
            ck = KSiOSi
            c0 = cos_SiOSi
         case default
         end select
         do kk = 1, ncmax(atom(j))
            k = proximity(j,kk)
            if (k == i .or. k == 0) cycle
            eang = energy_angle(rxyz(i,:),rxyz(j,:),rxyz(k,:),c0,ck)
            en = en + eang
            ! if (ANY(list(1:nl) == k)) cycle
            if (in_list(k,list,nl)) then
               !
               !  list   \
               !       i--\--j
               !           \/   j  outside list
               !           /\   BUT k inside list
               !          k  \
               ! if k is also in list then angle i-j-k will be counted twice
               ! so subtract off 1/2 the angle energy
               en = en - 0.5_wp*eang
            end if
         end do

      end do

      ! loop over all angles j-i-k
      !       k     centred around
      !      /      atom i
      ! j---i
      !
      select case(atom(i))
      case(iSilicon)
         ck = KOSiO
         c0 = cos_OSiO
      case(iOxygen)
         ck = KSiOSi
         c0 = cos_SiOSi
      case default
      end select
      do jj = 1,ncmax(atom(i))-1
         j = proximity(i,jj)
         if (j==0) cycle
         do kk = jj + 1,ncmax(atom(i))
            k = proximity(i,kk)
            if (k==0) cycle
            ! call angle_type(atom(j),atom(i),atom(k),c0,ck)
            en = en + energy_angle(rxyz(j,:),rxyz(i,:),rxyz(k,:),c0,ck)
         end do
      end do
      end do list_loop
   END FUNCTION keating_energy_List


!  PURE
   FUNCTION energy4(ifirst,ilast)
!-----function for the Keating energy
! for all atoms ifirst, ..., ilast
      USE coordinates_mod
      USE atom_types_mod, only: atom,ncmax
      USE connectivity_mod, only: proximity
      USE bond_angle_types_mod
      USE Keating_parameters_mod
      real(wp):: energy4
      integer,intent(in):: ifirst,ilast
      real(wp),parameter:: kbnd(2:4) = (/ KSiSi, KSiO, KSiO /)
      real(wp),parameter:: abnd(2:4) = (/ ASiSi, ASiO, ASiO /)
      real(wp),parameter:: acosang(1:3) = (/ cos_OSiO, cos_SiOSi, cos_SiOSi /)
      real(wp):: ebond,eangle,rij(NDIM),rik(NDIM),csa
      integer:: i,L,M,j,k
      real(wp):: kang(1:3,1:3,1:3)
      kang(2,1,2) = KOSiO
      kang(3,1,2) = KOSiO
      kang(2,1,3) = KOSiO
      kang(3,1,3) = KOSiO
      kang(1,2,1) = KSiOSi
      kang(1,3,1) = KSiOSi
      ebond = 0.0_wp
      eangle = 0.0_wp
!
      do i = ifirst,ilast

!     bond energy
      do L = 1,ncmax(atom(i))
         j = proximity(i,L)
         if (j == 0) cycle
         rij(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
         call pbc(rij)
         ebond = ebond + 0.25_wp*kbnd(atom(i)+atom(j)) &
               *(sqrt(dot_product(rij,rij)) - abnd(atom(i)+atom(j)))**2
      end do
!
!     angles with i as centre atom
      do L = 1,ncmax(atom(i))-1
         j = proximity(i,L)
         if (j == 0) cycle
         do M = L + 1,ncmax(atom(i))
            k = proximity(i,M)
            if (k == 0) cycle
            rij(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            rik(1:NDIM) = rxyz(k,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rij)
            call pbc(rik)
            csa = dot_product(rij,rik)/(sqrt(dot_product(rij,rij))*sqrt(dot_product(rik,rik)))
            eangle = eangle + 0.5_wp*kang(atom(j),atom(i),atom(k))*(csa - acosang(atom(i)))**2
         end do
      end do

      end do
      energy4 = ebond + eangle
   END FUNCTION energy4

END MODULE energy_keating_mod


!!>include 'force_keating_mod_sio2only.f90'

MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none

CONTAINS

   SUBROUTINE force_keating_bond(nbond,ibond,kbond,abond,en)
! Calculates the Keating Energy & Forces
! by summing over the bond & angle lists
      USE coordinates_mod
      integer,intent(in):: nbond,ibond(:,:)
      real(wp),intent(in):: kbond(:),abond(:)
      real(wp),intent(out):: en
      integer:: ib,ia1,ia2
      real(wp):: df(NDIM),kb,dist,r1(NDIM),r1sq
!
      en = 0.0_wp
      do ib = 1,nbond  ! Bond Stretch
         kb = kbond(ib)
         dist = Abond(ib)
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         r1sq = sqrt(dot_product(r1,r1))
         df = kb*(r1sq - dist)*r1/r1sq
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) - df
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) + df
         en = en + 0.5_wp*kb*(r1sq - dist)**2
      end do
   END SUBROUTINE force_keating_bond

   SUBROUTINE force_keating_angle(nang,iang,kang,cang,en)
! Calculates the total Keating Energy & Forces
! by summing over the bond & angle lists
      USE coordinates_mod
      integer,intent(in):: nang,iang(:,:)
      real(wp),intent(in):: kang(:),cang(:)
      real(wp),intent(out):: en
      integer:: ia,ia1,ia2,ia3
      real(wp):: r1(NDIM),r2(NDIM),df1(NDIM),df3(NDIM),r1sq,r2sq
      real(wp):: r1r2,ktheta,cos_theta,coef,coef1,coef2,coef3
!
      en = 0.0_wp
      do ia = 1,nang  ! Bond Angle
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         r2 = rxyz(ia3,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         call pbc(r2)
         r1sq = sqrt(dot_product(r1,r1))
         r2sq = sqrt(dot_product(r2,r2))
         r1r2 = dot_product(r1,r2)
         cos_theta = (r1r2)/(r1sq*r2sq)
         coef = ktheta*(cos_theta - cang(ia))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
         df1 = -coef*(r2/coef2 - r1r2*r1/coef1)
         df3 = -coef*(r1/coef2 - r1r2*r2/coef3)
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) + df1
         fxyz(ia3,1:NDIM) = fxyz(ia3,1:NDIM) + df3
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) - df1 - df3
         en = en + 0.5_wp*ktheta*(cos_theta - cang(ia))**2
      end do
   END SUBROUTINE force_keating_angle

   SUBROUTINE force_keating(en)
! Calculates the Keating Energy & Forces
! by summing over the bond & angle lists
      USE bond_angle_list_mod
      USE coordinates_mod
      real(wp),intent(out):: en
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: kb,dist,r1(NDIM),r1sq,ktheta,r2(NDIM),r2sq
      real(wp):: df(NDIM),df1(NDIM),df3(NDIM)
      real(wp):: r1r2,cos_theta,coef,coef1,coef2,coef3
!
      en = 0.0_wp
      do ib = 1,nbondtot  ! Bond Stretch
         kb = kbond(ib)
         dist = Abond(ib)
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         r1sq = sqrt(dot_product(r1,r1))
         df = kb*(r1sq - dist)*r1/r1sq
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) - df
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) + df
         en = en + 0.5_wp*kb*(r1sq - dist)**2
      end do
!
      do ia = 1,nang  ! Bond Angle
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         r1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         r2 = rxyz(ia3,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(r1)
         call pbc(r2)
         r1sq = sqrt(dot_product(r1,r1))
         r2sq = sqrt(dot_product(r2,r2))
         r1r2 = dot_product(r1,r2)
         cos_theta = (r1r2)/(r1sq*r2sq)
         coef = ktheta*(cos_theta - ctheta(ia))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
         df1 = -coef*(r2/coef2 - r1r2*r1/coef1)
         df3 = -coef*(r1/coef2 - r1r2*r2/coef3)
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) + df1
         fxyz(ia3,1:NDIM) = fxyz(ia3,1:NDIM) + df3
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) - df1 - df3
         en = en + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE force_keating


   PURE SUBROUTINE force_bond(r1,r2,d0,kb,en,df)
      USE coordinates_mod
      real(wp),intent(in):: r1(:),r2(:),d0,kb
      real(wp),intent(out):: en,df(NDIM)
      real(wp):: r12(NDIM),r12l
      r12 = r1 - r2
      call pbc(r12)
      r12l = sqrt(dot_product(r12,r12))
!     r12l = sqrt( r12(1)**2 + r12(2)**2 + r12(3)**2 )
      df = -kb*(r12l - d0)*r12/r12l
      en = 0.5_wp*kb*(r12l - d0)**2
   END SUBROUTINE force_bond


   PURE SUBROUTINE force_angle(r1,r2,r3,c0,ka,en,df1,df3)
      USE coordinates_mod
      real(wp),intent(in):: r1(:),r2(:),r3(:),c0,ka
      real(wp),intent(out):: en,df1(NDIM),df3(NDIM)
      real(wp):: r21(NDIM),r23(NDIM),cos_theta,r1sq,r2sq,r1r2
      real(wp):: coef,coef1,coef2,coef3
      r21 = r1 - r2
      r23 = r3 - r2
      call pbc(r21)
      call pbc(r23)
      r1sq = sqrt(dot_product(r21,r21))
      r2sq = sqrt(dot_product(r23,r23))
      r1r2 = dot_product(r21,r23)
      cos_theta = (r1r2)/(r1sq*r2sq)
      coef = ka*(cos_theta - c0)
      coef1 = r1sq**3*r2sq
      coef2 = r1sq*r2sq
      coef3 = r1sq*r2sq**3
      df1 = -coef*(r23/coef2 - r1r2*r21/coef1)
      df3 = -coef*(r21/coef2 - r1r2*r23/coef3)
      en = 0.5_wp*ka*(cos_theta - c0)**2
   END SUBROUTINE force_angle


!  PURE
   SUBROUTINE keating_force_atom(ia,en,df)
! Calculates the Keating energy of all the bonds & angles containing atom ia
! and the corresponding force on ia
      USE atom_types_mod
      USE coordinates_mod
      USE connectivity_mod
      USE bond_angle_types_mod
      USE Keating_parameters_mod
!     USE bond_angle_types_mod
      integer,intent(in):: ia
      real(wp),intent(out):: en,df(NDIM)
      integer:: j,jj,k,kk
      real(wp):: c0,ck !,kb,d0
      real(wp):: eang,ebnd,df1(NDIM),df3(NDIM)
!
      en = 0.0_wp
      df = 0.0_wp
      select case(atom(ia))
      case(iSilicon)
         ck = KOSiO
         c0 = cos_OSiO
      case(iOxygen)
         ck = KSiOSi
         c0 = cos_SiOSi
      case default
      end select
! the angles with ia as centre atom (j--ia--k)
      do jj = 1,ncmax(atom(ia))-1
         j = proximity(ia,jj)
         if (j == 0) cycle
         do kk = jj + 1,ncmax(atom(ia))
            k = proximity(ia,kk)
            if (k == 0) cycle
            ! call angle_type(atom(j),atom(ia),atom(k),c0,ck)
            call force_angle( rxyz(j,:),rxyz(ia,:),rxyz(k,:),c0,ck,eang,df1,df3 )
            df = df - df1 - df3
            en = en + eang
         end do
      end do

! bonds and other angles involving atom ia
      do jj = 1,ncmax(atom(ia))
         j = proximity(ia,jj)
         if ( j == 0) CYCLE
         ! call bond_type(atom(ia),atom(j),d0,kb)
         call force_bond( rxyz(ia,:),rxyz(j,:),ASiO,KSiO,ebnd,df1)
         df = df + df1
         en = en + ebnd
         select case(atom(j))
         case(iSilicon)
            ck = KOSiO
            c0 = cos_OSiO
         case(iOxygen)
            ck = KSiOSi
            c0 = cos_SiOSi
         case default
         end select
         do kk = 1,ncmax(atom(j))
            k = proximity(j,kk)
            if (k == 0) CYCLE
            if (k == ia) CYCLE
            ! call angle_type(atom(ia),atom(j),atom(k),c0,ck)
            call force_angle( rxyz(ia,:),rxyz(j,:),rxyz(k,:),c0,ck,eang,df1,df3 )
            df = df + df1
            en = en + eang
         end do
      end do
   END SUBROUTINE keating_force_atom


   SUBROUTINE keating_force_List(nl,list,en)
! Calculate the total bond & angle energy involving atoms 'i' in list
! interacting with atoms j & k either inside or outside the list
! and assign the force to the list atoms.
      USE atom_types_mod
      USE coordinates_mod
      USE connectivity_mod
      USE bond_angle_types_mod
      USE Keating_parameters_mod
      USE list_mod
      real(wp),intent(out):: en
      integer,intent(in):: nl,list(:)
      real(wp):: ck,c0 !,d0,kb
      real(wp):: eang,ebnd,df1(NDIM),df3(NDIM)
      integer:: i,il,j,jj,k,kk
!
      en = 0.0_wp
      list_loop: do il = 1, nl
      i = list(il)

      ! loop over atoms j bonded to i
      do jj = 1, ncmax(atom(i))

         j = proximity(i,jj)
         if (j == 0) cycle
         ! bond i-j
         ! call bond_type(atom(i),atom(j),d0,kb)
         call force_bond( rxyz(i,:),rxyz(j,:),ASiO,KSiO,ebnd,df1)
         fxyz(i,:) = fxyz(i,:) + df1
         en = en + 0.5_wp*ebnd

         ! if (ANY(list(1:nl) == j)) cycle
         if (in_list(j,list,nl)) cycle
         !
         ! if j is not in list
         !          |     k
         !  list    |   /
         !       i--|--j----k
         !          |   \
         !          |    k  outside list
         !
         ! then must include angles i-j-k
         ! and add the remaining 1/2 of the bond energy
         !
         en = en + 0.5_wp*ebnd
         ! call angle_type(atom(i),atom(j),atom(k),c0,ck)
         select case(atom(j))
         case(iSilicon)
            ck = KOSiO
            c0 = cos_OSiO
         case(iOxygen)
            ck = KSiOSi
            c0 = cos_SiOSi
         case default
         end select
         do kk = 1, ncmax(atom(j))
            k = proximity(j,kk)
            if (k == i .or. k == 0) cycle
            call force_angle( rxyz(i,:),rxyz(j,:),rxyz(k,:),c0,ck,eang,df1,df3 )
            fxyz(i,:) = fxyz(i,:) + df1
            en = en + eang
            ! if (ANY(list(1:nl) == k)) cycle
            if (in_list(k,list,nl)) then
               !
               !  list   \
               !       i--\--j
               !           \/   j  outside list
               !           /\   BUT k inside list
               !          k  \
               ! if k is also in list then angle i-j-k will be counted twice
               ! so subtract off 1/2 the angle energy
               en = en - 0.5_wp*eang
            end if
         end do

      end do

      ! loop over all angles j-i-k
      !       k     centred around
      !      /      atom i
      ! j---i
      !
      select case(atom(i))
      case(iSilicon)
         ck = KOSiO
         c0 = cos_OSiO
      case(iOxygen)
         ck = KSiOSi
         c0 = cos_SiOSi
      case default
      end select
      do jj = 1,ncmax(atom(i))-1
         j = proximity(i,jj)
         if (j==0) cycle
         do kk = jj + 1,ncmax(atom(i))
            k = proximity(i,kk)
            if (k==0) cycle
            ! call angle_type(atom(j),atom(i),atom(k),c0,ck)
            call force_angle( rxyz(j,:),rxyz(i,:),rxyz(k,:),c0,ck,eang,df1,df3 )
            fxyz(j,:) = fxyz(j,:) + df1
            fxyz(i,:) = fxyz(i,:) - df1 - df3
            fxyz(k,:) = fxyz(k,:) + df3
            en = en + eang
         end do
      end do
      end do list_loop
   END SUBROUTINE keating_force_List


   SUBROUTINE force4(ifirst,ilast,en)
!-----function for the Keating energy
! for all atoms ifirst, ..., ilast
      USE coordinates_mod
      USE atom_types_mod, only: atom,ncmax
      USE connectivity_mod, only: proximity
      USE Keating_parameters_mod
      USE bond_angle_types_mod
      real(wp),intent(out):: en
      integer,intent(in):: ifirst,ilast
      real(wp),parameter:: kbnd(2:4) = (/ KSiSi, KSiO, KSiO /)
      real(wp),parameter:: abnd(2:4) = (/ ASiSi, ASiO, ASiO /)
      real(wp),parameter:: acosang(1:3) = (/ cos_OSiO, cos_SiOSi, cos_SiOSi /)
      real(wp):: ebond,eangle,rij(NDIM),rik(NDIM),rijsq,kb,d0,c0,ck
      real(wp):: r1sq,r2sq,r1r2,cos_theta,coef,coef1,coef2,coef3,df(NDIM),df1(NDIM),df3(NDIM)
      integer:: i,L,M,j,k
      real(wp):: kang(1:3,1:3,1:3)
      kang(2,1,2)= KOSiO
      kang(3,1,2)= KOSiO
      kang(2,1,3)= KOSiO
      kang(3,1,3)= KOSiO
      kang(1,2,1)= KSiOSi
      kang(1,3,1)= KSiOSi
      ebond = 0.0_wp
      eangle = 0.0_wp
!
      do i = ifirst,ilast

!     bond energy
      do L = 1,ncmax(atom(i))
         j = proximity(i,L)
         if (j == 0) cycle
         kb = kbnd(atom(i)+atom(j))
         d0 = abnd(atom(i)+atom(j))
         rij(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
         call pbc(rij)
         rijsq = sqrt(dot_product(rij,rij))
         df = kb*(rijsq - d0)*rij/rijsq
         fxyz(i,1:NDIM) = fxyz(i,1:NDIM) + 0.5_wp*df
         fxyz(j,1:NDIM) = fxyz(j,1:NDIM) - 0.5_wp*df
         ebond = ebond + 0.25_wp*kb*(rijsq - d0)**2
!print '(2i6,2f12.6)',i,j,kbnd(atom(i)+atom(j)),abnd(atom(i)+atom(j))
      end do
!
!     angles with i as centre atom
      c0 = acosang(atom(i))
      do L = 1,ncmax(atom(i))-1
         j = proximity(i,L)
         if (j == 0) cycle
         do M = L + 1,ncmax(atom(i))
            k = proximity(i,M)
            if (k == 0) cycle
            ck = kang(atom(j),atom(i),atom(k))
            rij(1:NDIM) = rxyz(j,1:NDIM) - rxyz(i,1:NDIM)
            rik(1:NDIM) = rxyz(k,1:NDIM) - rxyz(i,1:NDIM)
            call pbc(rij)
            call pbc(rik)
            r1sq = sqrt(dot_product(rij,rij))
            r2sq = sqrt(dot_product(rik,rik))
            r1r2 = dot_product(rij,rik)
            cos_theta = (r1r2)/(r1sq*r2sq)
            coef = ck*(cos_theta - c0)
            coef1 = r1sq**3*r2sq
            coef2 = r1sq*r2sq
            coef3 = r1sq*r2sq**3
            df1 = -coef*(rik/coef2 - r1r2*rij/coef1)
            df3 = -coef*(rij/coef2 - r1r2*rik/coef3)
            fxyz(j,1:NDIM) = fxyz(j,1:NDIM) + df1
            fxyz(k,1:NDIM) = fxyz(k,1:NDIM) + df3
            fxyz(i,1:NDIM) = fxyz(i,1:NDIM) - df1 - df3
            eangle = eangle + 0.5_wp*ck*(cos_theta - c0)**2
         end do
      end do

      end do
      en = ebond + eangle
   END SUBROUTINE force4

END MODULE force_keating_mod

!!>include 'statistics_mod.f90'

MODULE statistics_mod
   USE precision_mod
   implicit none
!
! This module uses code modified from Alan Miller:
! (http://users.bigpond.net.au/amiller/update.f90)
! This is a collection of 3 routines to calculate the new mean
! and sum of squares of deviations about the mean (and hence the new
! sample variance or standard deviation) after:
!
! 1. A new observation becomes available (routine update).
! 2. An observation is dropped (routine downdate).
! 3. An observation is replaced with a new one (routine replace).
!
! Latest revision - 24 November 2001
! Alan Miller (amiller @ bigpond.net.au)
!
   TYPE statistic
      integer:: n=0
      real(wp):: mean=0.0_wp
      real(wp):: sumx=0.0_wp
      real(wp):: sumsq=0.0_wp
   END TYPE statistic

CONTAINS

   SUBROUTINE update(T, x)
!     x contains the value of the new case.
      type(statistic),intent(inout):: T
      real(wp),intent(in):: x
      real(wp):: dev
      T%n = T%n + 1
      T%sumx = T%sumx + x
      dev = x - T%mean
      T%mean = T%mean + dev/T%n
      T%sumsq = T%sumsq + dev*(x - T%mean)
   END SUBROUTINE update


   SUBROUTINE downdate(T, x)
!     x contained the value to be removed.
      type(statistic),intent(inout):: T
      real(wp),intent(in):: x
      real(wp):: dev
      T%n = T%n - 1
      T%sumx = T%sumx - x
      dev = x - T%mean
      T%mean = T%mean - dev/T%n
      T%sumsq = T%sumsq - dev*(x - T%mean)
   END SUBROUTINE downdate


   SUBROUTINE replace(T, x, y)
!     x is the value to be removed and replaced with the value y.
      type(statistic),intent(inout):: T
      real(wp),intent(in):: x, y
      real(wp):: diff, devx
      diff = y - x
      T%sumx = T%sumx + diff
      devx = x - T%mean
      T%mean = T%mean + diff/T%n
      T%sumsq = T%sumsq + diff*(devx + y - T%mean)
   END SUBROUTINE replace

   SUBROUTINE write_statistic(T,iu)
      integer,intent(in):: iu
      type(statistic),intent(in):: T
      write(iu,'(i9,3g18.8)') T%n, T%mean, T%sumx, T%sumsq
   END SUBROUTINE write_statistic

END MODULE statistics_mod

!!>include 'timer_mod.f90'

MODULE timer_mod
   USE statistics_mod
   implicit none

   TYPE timer
      real:: t_old
      real:: t_new
      type(statistic):: dt
   END TYPE timer

CONTAINS

   SUBROUTINE start_clock(T)
      type(timer):: T
      call cpu_time(T%t_old)
   END SUBROUTINE start_clock

   SUBROUTINE stop_clock(T)
      type(timer):: T
      real(wp):: delt
      call cpu_time(T%t_new)
      delt = T%t_new - T%t_old
      call update(T%dt, delt)
   END SUBROUTINE stop_clock

!   SUBROUTINE start_clock(T)
!      type(timer):: T
!      integer:: i0
!      call system_clock(i0)
!       T%t_old = real(i0,wp)
!   END SUBROUTINE start_clock
!
!   SUBROUTINE stop_clock(T)
!      type(timer):: T
!      real(wp):: delt
!      integer:: i1
!      call system_clock(i1)
!      T%t_new = real(i1,wp)
!      delt = T%t_new - T%t_old
!      call update(T%dt, delt)
!   END SUBROUTINE stop_clock

END MODULE timer_mod


!!>include 'force_energy_minimize_mod.f90'

MODULE FORCE_ENERGY_minimize_mod
   USE precision_mod
   implicit none
CONTAINS

   SUBROUTINE atom_energy(ia,en)
      USE energy_keating_mod
      USE repul_energy_mod
      integer,intent(in):: ia
      real(wp),intent(out):: en
      en = keating_energy_atom(ia) + repul_energy(ia,ia)
   END SUBROUTINE atom_energy

   SUBROUTINE force_solid(ifirst,ilast,En)
      USE force_keating_mod
      USE repul_energy_mod
      integer,intent(in):: ifirst,ilast
      real(wp),intent(out):: En
      real(wp):: ekeating, erepul
      call force4(ifirst,ilast,ekeating)
      call repul_force2(ifirst,ilast,erepul)
      en = ekeating + erepul
   END SUBROUTINE force_solid

   SUBROUTINE force_solid_list(nl,list,En)
      USE force_keating_mod
      USE repul_energy_mod
      integer,intent(in):: nl,list(:)
      real(wp),intent(out):: En
      real(wp):: ekeating, erepul
      call keating_force_List(nl,list,ekeating)
      call repul_force_List(nl,list,erepul)
      en = ekeating + erepul
   END SUBROUTINE force_solid_list

   SUBROUTINE Force1(ifirst,ilast,R,F,En)
      USE LJ_el_mod
      USE atom_types_mod
      integer,intent(in):: ifirst,ilast
      real(wp),intent(in):: R(:,:)
      real(wp),intent(inout):: F(:,:)
      real(wp),intent(out):: En
      call FORCE_ENERGY_LJ2(ifirst,ilast,atom,R,F,En)
   END SUBROUTINE Force1

   SUBROUTINE Force1_List(nl,list,R,F,En)
      USE LJ_el_mod
      USE atom_types_mod
      integer,intent(in):: nl,list(:)
      real(wp),intent(in):: R(:,:)
      real(wp),intent(inout):: F(:,:)
      real(wp),intent(out):: En
      call FORCE_ENERGY_LJ_LIST(nl,list,atom,R,F,En)
   END SUBROUTINE Force1_List

   SUBROUTINE ENERGY1(ia,ifirst,ilast,R,En)
      USE LJ_el_mod
      USE atom_types_mod
      integer,intent(in):: ia,ifirst,ilast
      real(wp),intent(in):: R(:,:)
      real(wp),intent(out):: En
      call ENERGY_LJ2(ia,ifirst,ilast,atom,R,En)
   END SUBROUTINE ENERGY1

   SUBROUTINE ENERGY1_List(ia,nl,list,R,En)
      USE LJ_el_mod
      USE atom_types_mod
      integer,intent(in):: ia,nl,list(:)
      real(wp),intent(in):: R(:,:)
      real(wp),intent(out):: En
      call ENERGY_LJ_LIST(ia,nl,list,atom,R,En)
   END SUBROUTINE ENERGY1_List

END MODULE FORCE_ENERGY_minimize_mod


!!>include 'Steepest_Descent_mod.f90'

MODULE  Steepest_Descent_mod
   USE precision_mod
   USE FORCE_ENERGY_minimize_mod, Force => force_solid, Force_List => force_solid_list
!  USE FORCE_ENERGY_minimize_mod, Force => Force1, Force_List=> Force1_List
   implicit none

CONTAINS

   SUBROUTINE Steepest_Descent(nstep,dRmax,ifirst,ilast,R,F)
      USE coordinates_mod
      USE nlist_mod
! USE BOX_FRAMES_MOD
      integer,intent(in):: nstep
      real(wp),intent(in):: dRmax
      integer,intent(in):: ifirst,ilast
      real(wp),intent(inout):: R(:,:),F(:,:)
      real(wp), parameter :: dev_rel = 1.0e-4_wp
!
! Minimize the energy
! Steepest Descent Algorithm
!
      real(wp),dimension(size(R,1),size(R,2)):: Ro,Fo
      real(wp):: Unew, Uold
      real(wp):: deltaR, delRF
      integer :: itry

      deltaR = dRmax
!     call Force(ifirst,ilast,R,F,Uold)
      call Force(ifirst,ilast,Uold)

! write (*,'(/a,1g18.8)') 'Initial Energy        : ',Uold

      do itry = 1, nstep

         Ro = R
         Fo = F

         delRF = deltaR/maxval(abs(F))
         call store_nlist

         R = R + delRF*F
         call pbc(R)
         call new_nlist(1,natom)

! call new_frame_file(imve,'frame',itry)
! call write_box_frame(imve,0,print_all = .true.)

!        call Force(ifirst,ilast,R,F,Unew)
         call Force(ifirst,ilast,Unew)

         if ( (abs((Unew - Uold)/Unew) < dev_rel) ) exit

         if (Unew < Uold) then
            Uold = Unew
            deltaR = deltaR*1.2_wp
            deltaR = min(deltaR , dRmax)
         else
            F = Fo
            R = Ro
            deltaR = deltaR*0.5_wp
            call restore_nlist
         end if

! write(*,'(a,i6,a,e12.4)') 'Step = ',itry,' :  Energy = ', Unew
! write(666,'(i6,4g18.8)') itry,Unew,Uold,delRF

      end do
! write (*,'(a,1g18.8/)') 'Final Energy          : ',Uold
   END SUBROUTINE Steepest_Descent


   SUBROUTINE Steepest_Descent_List(nstep,dRmax,nl,list,R,F)
      USE coordinates_mod
      USE nlist_mod
 USE BOX_FRAMES_MOD
      integer,intent(in):: nstep
      real(wp),intent(in):: dRmax
      integer,intent(in):: nl,list(:)
      real(wp),intent(inout):: R(:,:),F(:,:)
      real(wp), parameter :: dev_rel = 1.0e-4_wp
!
! Minimize the energy
! Steepest Descent Algorithm
!
      real(wp),dimension(nl,size(R,2)):: Ro,Fo
      real(wp):: Unew, Uold
      real(wp):: deltaR, delRF,maxabsf
      integer :: itry,il,ia

      deltaR = dRmax
!     call Force_List(nl,list,R,F,Uold)
      call Force_List(nl,list,Uold)

 write (*,'(/a,1g18.8)') 'Initial Energy        : ',Uold

      do itry = 1, nstep

         maxabsf = 0.0_wp
         do il = 1, nl
            ia = list(il)
            Ro(il,:) = R(ia,:)
            Fo(il,:) = F(ia,:)
            maxabsf = max(maxabsf,maxval(abs(F(ia,:))))
         end do
         delRF = deltaR/maxabsf
         ! delRF = deltaR/maxval(abs(F))
         call store_nlist

         do il = 1, nl
            ia = list(il)
            R(ia,:) = R(ia,:) + delRF*F(ia,:)
            call pbc(R(ia,:))
         end do
         call new_nlist(1,natom)

 call new_frame_file(imve,'frame',itry)
 call write_box_frame(imve,0,print_all = .true.)

!        call Force_List(nl,list,R,F,Unew)
         call Force_List(nl,list,Unew)

         if ( (abs((Unew - Uold)/Unew) < dev_rel) ) exit

         if (Unew < Uold) then
            Uold = Unew
            deltaR = deltaR*1.2_wp
            deltaR = min(deltaR , dRmax)
         else
            do il = 1, nl
               ia = list(il)
               F(ia,:) = Fo(il,:)
               R(ia,:) = Ro(il,:)
            end do
            deltaR = deltaR*0.5_wp
            call restore_nlist
         end if

 write(*,'(a,i6,a,e12.4)') 'Step = ',itry,' :  Energy = ', Unew
 write(666,'(i6,4g18.8)') itry,Unew,Uold,delRF

      end do
 write (*,'(a,1g18.8/)') 'Final Energy          : ',Uold
   END SUBROUTINE Steepest_Descent_List


   SUBROUTINE Steepest_Descent2(nstep,dRmax,ifirst,ilast,R,F)
      USE coordinates_mod
      USE nlist_mod
! USE BOX_FRAMES_MOD
      integer,intent(in):: nstep
      real(wp),intent(in):: dRmax
      integer,intent(in):: ifirst,ilast
      real(wp),intent(inout):: R(:,:),F(:,:)
      real(wp), parameter :: dev_rel = 1.0e-4_wp
!
! Minimize the energy
! Steepest Descent Algorithm
!
      real(wp),dimension(size(R,1),size(R,2)):: Ro,Fo
      real(wp):: Unew, Uold
      real(wp):: deltaR, delRF, RmsF, RmsF_o
      integer :: itry,nat

      nat = ilast - ifirst + 1
      deltaR = dRmax
!     call Force(ifirst,ilast,R,F,Uold)
      call Force(ifirst,ilast,Uold)
      RmsF = rms(F,nat)
      write (*,*) 'Initial Energy        : ',Uold

      do itry = 1, nstep

         RmsF_o = RmsF
         Ro = R
         Fo = F

         delRF = deltaR/RmsF

         call store_nlist
         R = R + delRF*F
         call pbc(R)
         call new_nlist(1,natom)

! call new_frame_file(imve,'frame',itry)
! call write_box_frame(imve,0,print_all = .true.)

!        call Force(ifirst,ilast,R,F,Unew)
         call Force(ifirst,ilast,Unew)

         if ( (abs((Unew - Uold)/Unew) < dev_rel) ) exit
         RmsF = rms(F,nat)

         if (Unew < Uold) then
            Uold = Unew
            deltaR = deltaR*1.2_wp
            deltaR = min(deltaR , dRmax)
         else
            F = Fo
            R = Ro
            RmsF = RmsF_o
            deltaR = deltaR*0.5_wp
            call restore_nlist
         end if

! write(*,'(a,i6,a,e12.4)') 'Step = ',itry,' :  Energy = ', Unew
! write(777,*) itry,Uold,delRF

      end do

! write (*,*) 'Final Energy          : ',Unew
   CONTAINS

   pure function rms(f,n)
      real(wp):: rms
      integer,intent(in):: n
      real(wp),intent(in):: f(:,:)
      real(wp):: sum2
      integer:: i
 !    rms = sqrt(sum( F*F ) )
 !    rms = sqrt(sum( F(1:n,1:NDIM)**2 ) )
      sum2 = 0.0_wp
      do i = 1 , n
         sum2 = sum2 + dot_product( F(i,:), F(i,:) )
      end do
      rms = sqrt(sum2/n)
   end function rms

   END SUBROUTINE Steepest_Descent2

END MODULE  Steepest_Descent_mod


!!>include 'relax_mc_mod.f90'

MODULE relax_mc_mod
   USE precision_mod, only: wp
   USE FORCE_ENERGY_minimize_mod !, CalcEnergy => atom_energy, CalcEnergy_List => Energy1_List
   implicit none
   public:: RELAX_MC, RELAX_MC_LIST

CONTAINS

   SUBROUTINE RELAX_MC(nr,ifirst,ilast,kbT)
      USE coordinates_mod
!     USE verlet_list_mod
      USE nlist_mod
      USE rand_mod
 USE BOX_FRAMES_MOD
      implicit none
      integer,intent(in):: nr,ifirst,ilast
      real(wp),intent(in):: kbT
      real(wp),save:: dr = 0.8_wp/Angstrom
      integer:: ir,ia,nfail,ntot !,ic,icn
      real(wp):: coprxyz(NDIM),Uold,Unew,dU_kT
      logical:: update
!
      nfail = 0
      ntot = 0
      DO ir = 1,NR
      DO ia = ifirst,ilast

!        if (update(ia)) call new_vlist

         coprxyz = rxyz(ia,:)

!        call CalcEnergy(ia,ifirst,ilast,Rxyz,Uold)
         call atom_energy(ia,Uold)

         rxyz(ia,1:NDIM) = coprxyz(1:NDIM) + (2.0_wp*randvec(NDIM) - 1.0_wp)*dr
         call pbc(rxyz(ia,:))

         update = .false.
         if (cell(rxyz(ia,:)) /= cell(coprxyz)) then
            call new_nlist(1,natom)
            update = .true.
         end if
!         if (update(ia)) call NEW_VLIST
!         if (update(ia)) then
!            write (6, *) 'ERROR: Displacement too large for Verlet '
!            STOP
!         end if
!
!        call CalcEnergy(ia,ifirst,ilast,Rxyz,Unew)
         call atom_energy(ia,Unew)
!
         dU_kT = (Unew - Uold)/kbT

         if (dU_kT <  0.0_wp) GOTO 22
         if (dU_kT > 50.0_wp) GOTO 11
         if (rand() < exp(-dU_kT)) GOTO 22
11       CONTINUE
            rxyz(ia,:) = coprxyz
            nfail = nfail + 1
            if (update) then
               call new_nlist(1,natom)
            end if
22       CONTINUE
         ntot = ntot + 1

      END DO
 call new_frame_file(imve,'frame',ir)
 call write_box_frame(imve,0,print_all = .true.)
 write(888, '(i8,4g18.8)') ir,Uold
      END DO
      if (real(nfail)/ntot < 0.5_wp) then
         dr = dr*1.02_wp
      else
         dr = dr/1.02_wp
      end if
!     write(999,'(f0.6,2(i7,1x),f0.6)') dr,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX_MC


   SUBROUTINE RELAX_MC_LIST(nr,nl,list,kbT)
      USE coordinates_mod
!     USE verlet_list_mod
      USE nlist_mod
      USE rand_mod
 USE BOX_FRAMES_MOD
      implicit none
      integer,intent(in):: nr,nl,list(:)
      real(wp),intent(in):: kbT
      real(wp),save:: dr = 0.8_wp/Angstrom
      integer:: ir,il,ia,nfail,ntot !,ic,icn
      real(wp):: coprxyz(NDIM),Uold,Unew,dU_kT
      logical:: update
!
      nfail = 0
      ntot = 0
      DO ir = 1,NR
      DO il = 1,NL

         ia = list(il)

!        if (update(ia)) call new_vlist

         coprxyz = rxyz(ia,:)

!        call CalcEnergy_List(ia,nl,list,Rxyz,Uold)
         call atom_energy(ia,Uold)

         rxyz(ia,1:NDIM) = coprxyz(1:NDIM) + (2.0_wp*randvec(NDIM) - 1.0_wp)*dr
         call pbc(rxyz(ia,:))

         update = .false.
         if (cell(rxyz(ia,:)) /= cell(coprxyz)) then
            call new_nlist(1,natom)
            update = .true.
         end if
!         if (update(ia)) call NEW_VLIST
!         if (update(ia)) then
!            write (6, *) 'ERROR: Displacement too large for Verlet '
!            STOP
!         end if
!
!        call CalcEnergy_List(ia,nl,list,Rxyz,Unew)
         call atom_energy(ia,Unew)
!
         dU_kT = (Unew - Uold)/kbT

         if (dU_kT <  0.0_wp) GOTO 22
         if (dU_kT > 50.0_wp) GOTO 11
         if (rand() < exp(-dU_kT)) GOTO 22
11       CONTINUE
            rxyz(ia,:) = coprxyz
            nfail = nfail + 1
            if (update) then
               call new_nlist(1,natom)
            end if
22       CONTINUE
         ntot = ntot + 1

      END DO
 call new_frame_file(imve,'frame',ir)
 call write_box_frame(imve,0,print_all = .true.)
 write(888, '(i8,4g18.8)') ir,Uold
      END DO
      if (real(nfail)/ntot < 0.5_wp) then
         dr = dr*1.02_wp
      else
         dr = dr/1.02_wp
      end if
!     write(999,'(f0.6,2(i7,1x),f0.6)') dr,ntot,nfail,real(nfail)/ntot
   END SUBROUTINE RELAX_MC_LIST

END MODULE relax_mc_mod




!!>include 'Henrys_law_calc_mod1.f90'

MODULE Henrys_law_calc_mod
    USE precision_mod
    implicit none

CONTAINS

   SUBROUTINE Henrys_law_calc(pr,rlower,rupper,ntrial,fover,rcutH,Uljel,Khenry,facc)
   !SUBROUTINE Henrys_law_calc(pr,rlower,rupper,ntrial,fover,Uljel,Khenry,facc)
       USE rand_mod
       USE nlist_mod
       USE quat2mat_mod
       USE ran_point_sphere
       USE global_vars_mod, only: kboltzT
       USE lj_el_mod
       USE probe_mol_mod
      type(probe_mol),intent(in):: pr
      real(wp),intent(in):: rlower(3),rupper(3)
      integer,intent(in):: ntrial
      real(wp),intent(in):: fover,rcutH
      real(wp),intent(out):: Uljel,Khenry,facc
!     real,parameter:: expon_max = log(huge(1.0_wp))*0.5_wp  ! not standard f95
      real(wp):: expon_max
      real::t0,t1,t2,t3
      type(probe_mol):: pr1
      real(wp):: aa(3,3),q(4),x,y,z,rnaccept,dUlj,dUel,dKHenry,expon
      integer:: it,nat,i
      logical:: rotate
!
!          ________(b)
!         /|      /|
!        / |     / |   Trial insertions in cube
!       /__|____/  |
!       |  |____|__|
!       | /     |  /
!       |/      | /
!       |_______|/     (a) is position rlower
!      (a)             (b) is position rupper
!
      expon_max = log(huge(1.0_wp))*0.5_wp
      nat = pr%n
      pr1 = pr
      rotate = ((size(pr%r,2) > 1) .and. nat > 1)
      rnaccept = 0
      Uljel = 0.0_wp
      KHenry = 0.0_wp
! open(unit=76,FILE= 'khenry.dat')
      insert_loop: do it = 1,ntrial
         x = rlower(1) + rand()*(rupper(1) - rlower(1))
         y = rlower(2) + rand()*(rupper(2) - rlower(2))
         z = rlower(3) + rand()*(rupper(3) - rlower(3))
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
!print*,'henrylaw calculations'
!print*,'i',i
!print*,pr1%r(1:3,i)
         end do
call cpu_time(t0)
         ! first check if there is a large overlap
         if (overlap(fover)) cycle insert_loop
         rnaccept = rnaccept + 1.0_wp
!print*,'rnaccept =',rnaccept          
call cpu_time(t1)
!print*,'overlap check',t1-t0
call cpu_time(t2) 
call ENERGY_LJ_EL(pr1,1,natom,dUlj,dUel)
!call ENERGY_LJEL_NLIST(pr1,rcutH,dUlj,dUel)
call cpu_time(t3)
!print*,'energycheck',t3-t2
         Uljel = Uljel + dUlj + dUel
         expon = -(dUlj + dUel)/kboltzT
         dKHenry = exp(min(expon,expon_max))
         KHenry = KHenry + dKHenry
!write(76,*)Uljel,expon,dKHenry,KHenry
      end do insert_loop
      Uljel = Uljel/real(ntrial,wp)
      KHenry = KHenry/real(ntrial,wp)
      facc = rnaccept/real(ntrial,wp)
!print*,'real(ntrial,wp)',real(ntrial,wp)
call delete_probe_mol(pr1)
      RETURN
   contains
      pure function overlap(fr)
         real(wp),intent(in):: fr
         logical:: overlap
         integer:: k,jj,j,ic,nc,ncell(125),neigh
         real(wp):: dr(3),sig2j
         overlap = .false.
         probe_atom_loop: do k = 1,nat
         ic = CELL( pr1%r(1:3,k) )
         call NEIGCELL(ic,1,neigh,ncell)
         
         cell_loop: do jj = 1,neigh

            nc = ncell(jj)
            if (nc == 0) cycle cell_loop
        
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               dr(1:3) = pr1%r(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               sig2j = sigLJ_2(atom(j))
               if (atom(j) == iSilicon) sig2j = 1.0_wp/Angstrom
               if (dot_product(dr,dr) < ( fr*( pr1%rad(k) + sig2j ) )**2) then
                  overlap = .true.
                  RETURN
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
         end do probe_atom_loop
      end function overlap
   END SUBROUTINE Henrys_law_calc


   SUBROUTINE Henrys_law_profile(probe,rlower,rupper,ntrial,fover,rcutH,dl,nb,grid)
!   SUBROUTINE Henrys_law_profile(probe,rlower,rupper,ntrial,fover,dl,nb,grid)
      USE probe_mol_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: rlower(3),rupper(3),dl,fover,rcutH
      integer,intent(in):: ntrial
      integer,intent(out):: nb(3)
      real(wp),intent(out):: grid(:,:,:)
      real(wp):: rl(3),ru(3),side(3),del(3),Khenry,facc,Uljel
!real::t0,t1
      integer:: i,j,k
      side = rupper-rlower
      nb = int(side/dl)
      del = side/nb
      do k = 1,nb(3)
         rl(3) = rlower(3) + (k-1)*del(3)
         ru(3) = rlower(3) + k*del(3)
         do j = 1,nb(2)
            rl(2) = rlower(2) + (j-1)*del(2)
            ru(2) = rlower(2) + j*del(2)
            do i = 1,nb(1)
               rl(1) = rlower(1) + (i-1)*del(1)
               ru(1) = rlower(1) + i*del(1)
!call cpu_time(t0)
!call Henrys_law_calc(probe,rl,ru,ntrial,fover,Uljel,Khenry,facc)
               call Henrys_law_calc(probe,rl,ru,ntrial,fover,rcutH,Uljel,Khenry,facc)
               grid(i,j,k) = Khenry
!call cpu_time(t1)
!print '(3i5,4g16.8)',i,j,k,Khenry,Uljel,facc,t1-t0
            end do
         end do
      end do
   END SUBROUTINE Henrys_law_profile


   SUBROUTINE Henrys_law_grid(probe,Accessible,grid,Nz,ntrial,fover,rcutH,output)
!  SUBROUTINE Henrys_law_grid(probe,Accessible,grid,ntrial,fover,output)
      USE grid_data_mod
      USE probe_mol_mod
      USE nlist_mod
      type(probe_mol),intent(in):: probe
      logical,intent(in):: Accessible(:,:,:)
      type(grid3D),intent(in):: grid
      integer,intent(in):: Nz
      integer,intent(in):: ntrial
      real(wp),intent(in):: fover,rcutH
      real(wp),intent(out):: output(:,:,:)
      real(wp):: rl(3),ru(3),facc,Uljel,Khenry ,t0,t1
  !    real::t0,t1
      integer:: ix,iy,iz,ic
      do iz = 1,Nz
         rl(3) = grid%rlower(3) + (iz-1)*grid%del(3)
         ru(3) = rl(3) + grid%del(3)
!print*,'insidehenrylawgrid'         
!print*,'iz=',iz
         do iy = 1,grid%nbin(2)
            rl(2) = grid%rlower(2) + (iy-1)*grid%del(2)
            ru(2) = rl(2) + grid%del(2)
            do ix = 1,grid%nbin(1)
               if (.NOT.Accessible(ix,iy,iz)) cycle
               rl(1) = grid%rlower(1) + (ix-1)*grid%del(1)
               ru(1) = rl(1) + grid%del(1)
 !check the cell number
! ic = icell(ix,iy,iz)
! print*,'ic= ',ic
 call cpu_time(t0)
               call Henrys_law_calc(probe,rl,ru,ntrial,fover,rcutH,Uljel,Khenry,facc)
call cpu_time(t1)
!print*, t1 -t0
!if (ic == 1)stop 'checck'
               output(ix,iy,iz) = Khenry
            end do
         end do
      end do
   END SUBROUTINE Henrys_law_grid

END MODULE Henrys_law_calc_mod

!!>include 'density_profile_mod.f90'

MODULE density_profile_mod
   USE precision_mod
   USE atom_types_mod
   USE grid_data_mod
   implicit none

   interface den_profile
      module procedure den_profile_grid, den_profile_arrs
   end interface den_profile

CONTAINS

   SUBROUTINE init_density_profile(rlower,rupper,dl,ntypes,nden,ndenz)
      real(wp),intent(in):: rlower(3),rupper(3),dl
      integer,intent(in):: ntypes
      real(wp),allocatable,intent(inout):: nden(:,:,:,:),ndenz(:,:)
      real(wp):: side(3)
      integer:: nbin(3)
      side = rupper-rlower
      nbin = int(side/dl)
      allocate( nden(nbin(1),nbin(2),nbin(3),ntypes), ndenz(nbin(3),ntypes) )
   END SUBROUTINE init_density_profile


   SUBROUTINE den_profile_z0(ibegin,iend,Zlower,Zupper,dl,ndenz)
      USE coordinates_mod
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: Zlower,Zupper,dl
      real(wp),intent(out):: ndenz(:,:)
      real(wp):: side,delz,delzi
      integer:: j,iz,nbz
      side = Zupper-Zlower
      nbz = int(side/dl)
      delz = side/nbz
      delzi = 1.0_wp/delz
      ndenz = 0.0_wp
      do j = ibegin,iend
         iz = int((rxyz(j,3)-Zlower)*delzi) + 1
         if ((iz < 1).or.(iz > nbz)) cycle
         ndenz(iz,atom(j)) = ndenz(iz,atom(j)) + 1.0_wp
      end do
   END SUBROUTINE den_profile_z0


   SUBROUTINE write_calc_den_profile_z(iu,nstep,ibegin,iend,dl)
      use atom_types_mod
      use coordinates_mod
      integer,intent(in):: iu,nstep,ibegin,iend
      real(wp),intent(in):: dl
      real(wp),allocatable:: nden(:,:)
      real(wp):: delz,delzi,zvol,zmin,zmax,zj,zcline
      integer:: i,j,ii,nbinm,nbinp
!
! Calculate and print the z density profile in units of atoms/A^3
!
      zmax = maxval(rxyz(ibegin:iend,3))
      zmin = minval(rxyz(ibegin:iend,3))
      nbinm = int((zmin)/dl) - 2
      nbinp = int((zmax)/dl) + 2
      allocate( nden(nbinm:nbinp,7) )
!
      delz = dl
      delzi = 1.d0/delz
      nden = 0.0
      do j = ibegin,iend
         zj = rxyz(j,3)
         if (zj >= 0.0) then
            ii = int(zj*delzi) + 1
         else
            ii = int(zj*delzi) - 1
         end if

         if(atom(j) == iSilicon) then
            nden(ii,1) = nden(ii,1) + 1.0_wp
         else if (atom(j) == iOxygen) then
            nden(ii,2) = nden(ii,2) + 1.0_wp
         else if (atom(j) == iOxygenH) then
            nden(ii,3) = nden(ii,3) + 1.0_wp
         else if (atom(j) == iOxygenC) then
            nden(ii,4) = nden(ii,4) + 1.0_wp
         else if (atom(j) == iCarbonH2) then
            nden(ii,5) = nden(ii,5) + 1.0_wp
         else if (atom(j) == iCarbonH3) then
            nden(ii,6) = nden(ii,6) + 1.0_wp   
         end if
      end do
!
      zvol = boxl(1)*boxl(2)*delz*(Angstrom**3)
      write(iu,'(a,i6)') '#nstep ',nstep
      do i = nbinm,-1
         zcline = (i + 0.5_wp)*delz
         write(iu,'(7(2x,f0.8))') zcline,nden(i,1)/zvol,nden(i,2)/zvol,nden(i,3)/zvol, &
         nden(i,4)/zvol,nden(i,5)/zvol,nden(i,6)/zvol
      end do
      do i = 1,nbinp
         zcline = (i-0.5_wp)*delz
         write(iu,'(7(2x,f0.8))') zcline,nden(i,1)/zvol,nden(i,2)/zvol,nden(i,3)/zvol, &
         nden(i,4)/zvol,nden(i,5)/zvol,nden(i,6)/zvol
      end do
      write(iu,'(/)')
      call flush(iu)
      deallocate(nden)
   END SUBROUTINE write_calc_den_profile_z


   SUBROUTINE den_profile_arrs(ibegin,iend,rlower,rupper,dl,nden,ndenz)
      USE coordinates_mod
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: rlower(3),rupper(3),dl
      real(wp),intent(out):: nden(:,:,:,:),ndenz(:,:)
      real(wp):: side(3),delx,dely,delz
      real(wp):: delxi,delyi,delzi
      integer:: j,jj,ix,iy,iz,nbx,nby,nbz
      side = rupper - rlower
      nbx = int(side(1)/dl)
      nby = int(side(2)/dl)
      nbz = int(side(3)/dl)
      delx = side(1)/nbx
      dely = side(2)/nby
      delz = side(3)/nbz
      delxi = 1.0_wp/delx
      delyi = 1.0_wp/dely
      delzi = 1.0_wp/delz
      nden = 0.0_wp
      ndenz = 0.0_wp
      do j = ibegin,iend
         ix = int((rxyz(j,1)-rlower(1))*delxi) + 1
         iy = int((rxyz(j,2)-rlower(2))*delyi) + 1
         iz = int((rxyz(j,3)-rlower(3))*delzi) + 1
         if ((ix < 1).or.(ix > nbx)) cycle
         if ((iy < 1).or.(iy > nby)) cycle
         if ((iz < 1).or.(iz > nbz)) cycle
         jj = atom(j)
         nden(ix,iy,iz,jj) = nden(ix,iy,iz,jj) + 1.0_wp
      end do
      forall(iz=1:nbz,jj=1:ntyp) ndenz(iz,jj) = sum(nden(:,:,iz,jj))
   END SUBROUTINE den_profile_arrs


   SUBROUTINE den_profile_grid(ibegin,iend,G,nden,ndenz)
      USE coordinates_mod
      USE grid_data_mod
      integer,intent(in):: ibegin,iend
      type(grid3d),intent(in):: G
      real(wp),intent(out):: nden(:,:,:,:),ndenz(:,:)
      integer:: i,j,jj,ib(3)
      nden = 0.0_wp
      do j = ibegin,iend
         ib = int((rxyz(j,1:NDIM)-G%rlower)*G%deli) + 1
         if ((any(ib < 1).or.any(ib > G%nbin))) cycle
         jj = atom(j)
         nden(ib(1),ib(2),ib(3),jj) = nden(ib(1),ib(2),ib(3),jj) + 1.0_wp
      end do
!      nden = nden/G%volbin
      forall(i=1:G%nbin(3),jj=1:ntyp) ndenz(i,jj) = sum(nden(:,:,i,jj))   !/(G%nbin(1)*G%nbin(2))
   END SUBROUTINE den_profile_grid


   SUBROUTINE write_den_profile_z(iu,nstep,zlower,delz,ndenz)
      USE coordinates_mod
      integer,intent(in):: iu,nstep
      real(wp),intent(in):: zlower,delz,ndenz(:,:)
      integer:: i,j,itop
      do itop = size(ndenz,1),1,-1
         if (any(ndenz(itop,:) > 0.0_wp)) exit
      end do
      write(iu,'(a,i6)') '#nstep ',nstep
      write(iu,'(a)') '#density [atoms/A^3]'
      do i = 1,itop
         write(iu,'(20f12.6)') zlower+(i-0.5_wp)*delz,(ndenz(i,j)/(boxl(1)*boxl(2)*delz*angstrom**3),j=1,size(ndenz,2))
      end do
      write(iu,'(/)')
   END SUBROUTINE write_den_profile_z

END MODULE density_profile_mod

!!>include 'attach_molecule_mod.f90'

MODULE attach_molecule_mod
    USE precision_mod
    USE global_vars_mod
    USE constants_mod
    USE atom_types_mod
    USE coordinates_mod
    USE connectivity_mod
    USE rand_mod
    implicit none

CONTAINS

   SUBROUTINE attach_mol(L,ndel,idel,bondl,bondang,natt,rmol,conect,atype,success)
!  SUBROUTINE attach_mol(L,bondl,bondang,pmol,success)
      USE nlist_mod
      USE rotate_axis_mod
      USE verlet_list_mod
      USE HKNonLattice_mod
      integer,intent(inout):: L
      integer,intent(in):: ndel,idel(:)
      real(wp),intent(in):: bondang,bondl
!     type(probe_mol):: pmol
      integer,intent(in):: natt
      real(wp),intent(in):: rmol(:,:)
      integer,intent(in):: conect(:,:)
      integer,intent(in):: atype(:)
      logical,intent(out):: success
!
      real(wp):: xx,yy,zz,phi,cphi,sphi,cosang,sinang
      real(wp):: r2
      real(wp):: r3v(3),r12(3),r23(3),rSi(3),aa(3,3)
      integer:: i,j,LSi,ii,jj,M,ic,nc,nl,natom_old
      integer:: copproximity(natom_max,size(proximity,dim=2)),copatom(natom_max)
      real(wp):: coprxyz(natom_max,3)
      integer:: neighbors(343),nn
      logical:: overlap
!
      cosang = cos(bondang)
      sinang = sin(bondang)
      success = .false.
! save all the old data before the attempted attachment
      call NEW_VLIST(1,natom)
      do i = 1,natom
         coprxyz(i,1:3) = rxyz(i,1:3)
         copproximity(i,:) = proximity(i,:)
         copatom(i) = atom(i)
         cnlist(i) = nlist(i)
         clist(i,1:nlist(i)) = list(i,1:nlist(i))
      end do
      natom_old = natom
!
!-----Attaching a new molecule
      do j = 1,ncmax(atom(L))
         ii = proximity(L,j)
         if ( atom(ii) == iSilicon ) then
            LSi = ii
            exit
         end if
      end do
      if ( j > ncmax(atom(L)) ) stop 'error: j > ncmax(atom(L)'
      r3v(1:3) = rxyz(L,1:3) - rxyz(LSi,1:3)
      call pbc(r3v)
      r12 = r3v/sqrt(dot_product(r3v,r3v))

      do i = 1,ndel
         call delete_atom(idel(i))
      end do
      if ( L > natom) L = natom

!     prepare the position for silica
      call zaxis2vect(r12, aa)
      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      zz = -bondl*cosang
      xx = bondl*sinang*cphi
      yy = bondl*sinang*sphi
      r23 = matmul(aa,(/xx,yy,zz/))
      rSi(1:3) = rxyz(L,1:3) + r23
      r23 = r23/sqrt(dot_product(r23,r23))

!     prepare the position of the molecule
      call zaxis2vect(r23,aa)
!
!     attach the molecule(mol)
      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      do i = 1,natt
         rxyz(natom + i,1) = cphi*rmol(i,1) - sphi*rmol(i,2)
         rxyz(natom + i,2) = sphi*rmol(i,1) + cphi*rmol(i,2)
         rxyz(natom + i,3) = rmol(i,3)
      end do
      do i = 1,natt
         rxyz(natom + i,1:3) = matmul(aa,rxyz(natom + i,1:3))
      end do
      do i = 1,natt
         rxyz(natom + i,1) = rxyz(natom + i,1) + rSi(1)
         rxyz(natom + i,2) = rxyz(natom + i,2) + rSi(2)
         rxyz(natom + i,3) = rxyz(natom + i,3) + rSi(3)
      end do
!
!-----Apply PBC to the new atoms
      call pbc(rxyz(natom + 1:natom + natt,:))
!
      proximity(natom + 1:natom + natt,:) = conect(1:natt,:)
      where (proximity(natom + 1:natom + natt,:) /= 0) &
             proximity(natom + 1:natom + natt,:) = &
             proximity(natom + 1:natom + natt,:) + natom
      where (proximity(natom + 1,:) > (natom + natt)) proximity(natom + 1,:) = 0
      atom(natom + 1:natom + natt) = atype(1:natt)

      call set_proximity(natom + 1,0,L)
      call set_proximity(L,0,natom + 1)

      atom(L) = iOxygen
      atom(natom + 1:natom + natt) = atype(1:natt)
!
      atomL(natom + 1:natom + natt) = atomL(L)  ! cluster label for new atoms
      natom = natom + natt
      call NEW_NLIST(1,natom)
!
! Check for overlap
      overlap = .false.
      outer_atom_loop: do ii = natom - natt + 1,natom
      nl = nlayers_ll(sigLJ(atom(ii)))
      if (nl == 0) cycle
      ic = CELL(rxyz(ii,:))  !link list
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
         nc = neighbors(jj)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if ( M > natom - natt) GOTO 100
            if (atom(M) == iSilicon) GOTO 100
            if (M == ii) GOTO 100
            if (nearest_neighbor2(ii,m)) GOTO 100
            r3v(1:3) = rxyz(M,1:3) - rxyz(ii,1:3)
            call pbc(r3v)
            r2 = r3v(1)**2 + r3v(2)**2 + r3v(3)**2
            if (r2 <= (sigLJ_2(atom(m)) + sigLJ_2(atom(ii)))**2) then
               overlap = .true.
               exit outer_atom_loop
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
!
      if (.not.overlap) then
         success = .true.
         RETURN
      end if
!
999   CONTINUE
      atom(natom + 1:natom + natt) = 0
      proximity(natom + 1:natom + natt,:) = 0
      natom = natom_old
      do i = 1,natom
         rxyz(i,1:3) = coprxyz(i,1:3)
         proximity(i,:) = copproximity(i,:)
         atom(i) = copatom(i)
         nlist(i) = cnlist(i)
         list(i,1:cnlist(i)) = clist(i,1:cnlist(i))
      end do
      call NEW_NLIST(1,natom)
      RETURN
   END SUBROUTINE attach_mol

END MODULE attach_molecule_mod


!!>include 'voidage_mod.f90'

MODULE voidage_mod
   USE precision_mod
   implicit none
CONTAINS

   SUBROUTINE voidage_calc(pr,rlower,rupper,ntrial,voidage)
      USE coordinates_mod
      USE atom_types_mod
      USE constants_mod
      USE rand_mod
      USE nlist_mod
      USE quat2mat_mod
      USE ran_point_sphere
      USE probe_mol_mod
      type(probe_mol),intent(in):: pr
      real(wp),intent(in):: rlower(3),rupper(3)
      integer,intent(in):: ntrial
      real(wp),intent(out):: voidage
      real(wp),allocatable:: probe(:,:)
      real(wp):: aa(3,3),q(4)
      real(wp):: x,y,z,dr(3),rnaccept,sig2j
      integer:: i,k,j,ic,jj,nc,it,nat,neigh,ncell(125)
      logical:: rotate
      nat = pr%n
      allocate(probe(3,nat))
      rotate = ((size(probe,2) > 1))
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         x = rlower(1) + rand()*(rupper(1) - rlower(1))
         y = rlower(2) + rand()*(rupper(2) - rlower(2))
         z = rlower(3) + rand()*(rupper(3) - rlower(3))
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
         do i = 1,nat
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
                  dr(1:3) = probe(1:3,k) - rxyz(j,1:3)
                  call pbc(dr)
                  sig2j = sigLJ_2(atom(j))
                  if (sig2j < 0.01_wp/angstrom) sig2j = 0.6_wp/angstrom
                  if (dot_product(dr,dr) < (pr%rad(k) + sig2j)**2) then
                     cycle insert_loop
                  end if
                  j = LL(j)
               end do cell_atom_loop
            end do cell_loop
         end do probe_atom_loop
         rnaccept = rnaccept + 1.0_wp
      end do insert_loop
      voidage = rnaccept/real(ntrial,wp)
   END SUBROUTINE voidage_calc


   SUBROUTINE voidage_calc1(rprobe,rlower,rupper,ntrial,voidage)
      USE coordinates_mod
      USE atom_types_mod
      USE constants_mod
      USE rand_mod
      USE nlist_mod
      real(wp),intent(in):: rprobe,rlower(3),rupper(3)
      integer,intent(in):: ntrial
      real(wp),intent(out):: voidage
      real(wp):: r(3),dr(3),rnaccept,sig2j
      integer:: j,ic,jj,nc,it,neigh,ncell(125)
      rnaccept = 0
      insert_loop: do it = 1,ntrial
         r(1) = rlower(1) + rand()*(rupper(1) - rlower(1))
         r(2) = rlower(2) + rand()*(rupper(2) - rlower(2))
         r(3) = rlower(3) + rand()*(rupper(3) - rlower(3))
         ic = CELL(r)
         call NEIGCELL(ic,1,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0)cycle cell_loop
            j = HOC(nc)
            do while (j /= 0)
               dr(1:3) = r(1:3) - rxyz(j,1:3)
               call pbc(dr)
               sig2j = sigLJ_2(atom(j))
               if (sig2j < 0.01_wp/angstrom) sig2j = 0.6_wp/angstrom
               if (dot_product(dr,dr) < (rprobe + sig2j)**2) then
                  cycle insert_loop
               end if
               j = LL(j)
            end do
         end do cell_loop
         rnaccept = rnaccept + 1.0_wp
      end do insert_loop
      voidage = rnaccept/real(ntrial,wp)
   END SUBROUTINE voidage_calc1


   SUBROUTINE voidage_profile(probe,rlower,rupper,ntrial,dl,nb,voidage)
      USE coordinates_mod
      USE probe_mol_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: rlower(3),rupper(3),dl
      integer,intent(in):: ntrial
      integer,intent(out):: nb(3)
      real(wp),intent(out):: voidage(:,:,:)
      real(wp):: void
      real(wp):: rl(3),ru(3),side(3)
      real(wp):: del(3)
!     real::t0,t1
      integer:: i,j,k
      side = rupper-rlower
      nb = int(side/dl)
      del = side/nb
      do k = 1,nb(3)
         rl(3) = rlower(3) + (k-1)*del(3)
         ru(3) = rlower(3) + k*del(3)
         do j = 1,nb(2)
            rl(2) = rlower(2) + (j-1)*del(2)
            ru(2) = rlower(2) + j*del(2)
            do i = 1,nb(1)
               rl(1) = rlower(1) + (i-1)*del(1)
               ru(1) = rlower(1) + i*del(1)
               !call cpu_time(t0)
               call voidage_calc(probe,rl,ru,ntrial,void)
               voidage(i,j,k) = void
               !call cpu_time(t1)
               !print '(3i5,4g16.8)',k,j,i,void,t1-t0
            end do
         end do
      end do
   END SUBROUTINE voidage_profile

END MODULE voidage_mod




PROGRAM SIMBOX
      USE precision_mod
      USE phys_const_mod
      USE commandline_mod
      USE command_arg_mod
      USE math_const_mod
      USE phys_const_mod
      USE precision_mod
      USE utility_mod
      USE files_mod
      USE constants_mod
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE bond_angle_list_mod
      USE rand_mod
      USE atom_types_mod
      USE list_mod
      USE listset_mod
      USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE nlist_mod
      USE verlet_list_mod
      USE repul_energy_mod
      USE charges_mod
      USE force_en_el_mod
      USE insert_mod
      USE tcf_mod
      USE msd_mod
      USE vcf_mod
      USE rdf_mod
      USE Steepest_Descent_mod
      USE relax_mc_mod
      USE energy_keating_mod
      USE statistics_mod
      USE timer_mod
      USE find_atoms_mod
      USE probe_mol_mod
      USE probe_mol_init_mod
      USE check_structure_mod
      USE voidage_mod
      USE molecule_type_mod
      USE charmm_mod
      USE Dreiding_parameters_mod
      USE grid_data_mod
      USE HKNonLattice_mod
      USE hk_lattice_uf_mod
      USE PERCOLATION_LATTICE_MOD
      USE Henrys_law_calc_mod
      USE rotate_axis_mod
      USE tetra_coords_mod
      USE density_profile_mod

      implicit none
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      integer:: number_moves, nequib, nsamp, nprint_all, nfreq_movie, n_short
      logical:: make_movie
      integer:: ntau_msd, ntau_vcf
      real(wp):: del_nlist, fover,rl(NDIM),ru(NDIM),energy
      real(wp):: dt,qj,sumqj,T_kelvin,T_kinetic,sumx,sumy,sumz
      real(wp),allocatable:: mol_vec(:,:),mol_vec_h2o(:,:)
      real(wp):: EK,Ek_gas,Ek_O2,Ek_N2,Ek_CO2,Ek_tot,Ek_h2o
      real(wp):: dtreal,time,rr(NDIM),mtmp,ULJAtt,Ukeating,Ubond,Uang
      real(wp),allocatable:: ULJEL(:)
      real(wp),allocatable::e_bond_CO2(:),e_angle_CO2(:),e_bond_gas(:),e_angle_gas(:)
      real(wp):: e_bond_O2,e_bond_N2,Uold,Unew
      real(wp),allocatable:: massi(:,:),mass(:,:)
      real(wp),allocatable:: h2o_massi(:,:,:),h2o_mass(:,:,:),rad_h2o(:)
      real(wp),allocatable:: GAS_massi(:,:,:),GAS_mass(:,:,:),rad_GAS(:)
      real(wp),allocatable:: Ox_massi(:,:,:),Ox_mass(:,:,:)
      real(wp),allocatable:: N2_massi(:,:,:),N2_mass(:,:,:)
      real(wp),allocatable:: CO2_massi(:,:,:),CO2_mass(:,:,:)
      real(wp):: rad_Ox(3),rad_N2(3),rad_CO2(3)
      integer,allocatable:: mol_atoms(:), mol_atoms_h2o(:)
      integer,allocatable:: listOH(:),listALL(:),listSph(:),alst(:)
      integer:: nSph,nlst,n_oh,n_si,n_oc,n_ch3,n_ch2,n_ob,ntot !,nfixed
      integer:: j,i,istep,stat,nc
      integer:: itmp,narg,Ndf
      real(wp),allocatable:: coprxyz(:,:)
      real:: tp(0:10),sumtp(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5,c32
      character(len=132):: coordfile,restart_coordfile,restart_velfile,carg,ftyp*3,ctrlfile
      integer:: io_temp, io_energy, io_config, io_config_final, io_topol, io_ctrl
      integer:: io_restart_xyz, io_restart_v
      logical:: restart,ok,check_self_overlap=.false.
      type(tcf_type):: msd_gas, vcf_gas, msd_o2, vcf_o2, msd_n2, vcf_n2,msd_co2, vcf_co2, msd_h2o, vcf_h2o
      type(timer):: t(30)
      real(wp):: dLflexible, dLfixed, dlkeep,voidage(5),volbox,volbox1,dens_gcm3,fdens_wt
      integer:: nbinv(3)
      real(wp)::Rbox_min(3),Rbox_max(3),Rbox_min1(3),Rbox_max1(3),boxz,zmax,zmin,zmax1,zmin1
      real(wp),allocatable:: henryk(:,:,:,:)
      logical,allocatable:: Accessible(:,:,:,:)
      integer,allocatable:: Lattice(:,:,:,:)
      integer,parameter:: n_vapor_species = 1
      integer,save:: iCO2 = 1
      integer:: ncluster(n_vapor_species),n_try_Henry      
      type(mylist)::NumAcc(n_vapor_species)
      integer:: nbinmax,nz_top_atom,nz_top
      real(wp):: top_atom,void_crit,Rbox_top(3),fover_henry,del_bin,z_top,rcutH
      Integer:: nx,ny,nz,volbin
      type(grid3D):: Grid
      type(mylist):: fixed_list
      real(wp),allocatable:: nden(:,:,:,:),ndenz(:,:)
      
!
      print *,'Molecular Dynamics'
      print *,'Number of dimensions = ',NDIM
      print *,'Number of periodic dimensions = ',NPER
      if (NDIM < NPER) then
         print *,'The number of dimensions must be >= the number of periodic dimensions'
         STOP 'Change NDIM &/ NPER in coordinates_mod and recompile'
      end if
      call cpu_time(tp(0))
      call get_free_file_unit(io_energy)
      open(unit=io_energy,file='energy.out')
      call get_free_file_unit(io_temp)
      open(unit=io_temp,file='temperature.out')
!
      narg = command_argument_count()
      call next_command_argument(carg, "command",stat)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) trim(carg),'  control_file  atom_file  T[K]  {xyz_restart  vel_restart} '
         stop
      end if
!
      call next_command_argument(ctrlfile, "control file",stat)
      call get_free_file_unit(io_ctrl)
      open(unit=io_ctrl,file=trim(ctrlfile))
      read(io_ctrl,*) c32, number_moves        ; write(*,*) trim(c32), number_moves
      read(io_ctrl,*) c32, n_short             ; write(*,*) trim(c32), n_short
      read(io_ctrl,*) c32, nequib              ; write(*,*) trim(c32), nequib
      read(io_ctrl,*) c32, dtreal              ; write(*,*) trim(c32), dtreal
      read(io_ctrl,*) c32, nsamp               ; write(*,*) trim(c32), nsamp
      read(io_ctrl,*) c32, nprint_all          ; write(*,*) trim(c32), nprint_all
      read(io_ctrl,*) c32, make_movie          ; write(*,*) trim(c32), make_movie
      read(io_ctrl,*) c32, nfreq_movie         ; write(*,*) trim(c32), nfreq_movie
      read(io_ctrl,*) c32, ntau_msd            ; write(*,*) trim(c32), ntau_msd
      read(io_ctrl,*) c32, ntau_vcf            ; write(*,*) trim(c32), ntau_vcf
      read(io_ctrl,*) c32, n_O2,n_atom_o2      ; write(*,*) trim(c32), n_O2,n_atom_o2
      read(io_ctrl,*) c32, n_N2,n_atom_n2      ; write(*,*) trim(c32), n_N2,n_atom_n2
      read(io_ctrl,*) c32, n_CO2,n_atom_CO2    ; write(*,*) trim(c32), n_CO2,n_atom_CO2
      read(io_ctrl,*) c32, n_GAS, n_atom_gas   ; write(*,*) trim(c32), n_GAS, n_atom_gas
      read(io_ctrl,*) c32, n_h2o, n_atom_h2o   ; write(*,*) trim(c32), n_h2o, n_atom_h2o
      read(io_ctrl,*) c32, del_nlist           ; write(*,*) trim(c32), del_nlist
      read(io_ctrl,*) c32, fover               ; write(*,*) trim(c32), fover
      read(io_ctrl,*) c32, n_try_Henry         ; write(*,*) trim(c32), n_try_Henry
      read(io_ctrl,*) c32, fover_henry         ; write(*,*) trim(c32), fover_henry
      read(io_ctrl,*) c32,void_crit            ; write(*,*) trim(c32), void_crit

!      
      del_nlist = del_nlist/Angstrom
      del_bin = del_nlist
      rcutH = 10.0_wp/Angstrom ! not used in this prog

      call next_command_argument(coordfile, "input file",stat)
      nc = len_trim(coordfile)
print *,"'",trim(coordfile),"'"
print *,nc
      ftyp = coordfile(nc - 2:nc)
      call get_free_file_unit(io_config)
      open(unit=io_config,file=trim(coordfile))
      call next_command_argument(T_kelvin, "Temperature [K]",stat)
      kBoltzT = K_ev*T_kelvin
      call next_command_argument(restart_coordfile, "restart coord file",stat)
      if (stat /= 0) then
         restart = .false.
      else
         restart = .true.
         write (*,*) 'restarting from ', trim(restart_coordfile)
         call get_free_file_unit(io_restart_xyz)
         open(unit=io_restart_xyz,file=trim(restart_coordfile))
         call next_command_argument(restart_velfile, "restart vel file", stat)
         if (stat /= 0) then
            print *,'restart velocity file missing'
            stop
         else
            call get_free_file_unit(io_restart_v)
            open(unit=io_restart_v,file=trim(restart_velfile))
         end if
      end if

!
    
      dt = dtreal*dt_fs
print*,'))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))'
      call init_probe_mols()
print*,'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
      rad_Ox = (/ sig_O_O2,sig_O_O2,0.0_wp /)*0.5_wp
      rad_N2 = (/ sig_n_N2,sig_n_N2,0.0_wp /)*0.5_wp
      rad_CO2 = (/ sig_C_CO2,sig_O_CO2,sig_O_CO2 /)*0.5_wp

      allocate(mol_atoms(n_atom_gas))
      allocate(GAS_atom(1:n_GAS,n_atom_gas))
      allocate(GAS_charge(n_GAS,n_atom_gas))
      allocate(GAS_mass(1:n_GAS,n_atom_gas,3),GAS_massi(1:n_GAS,n_atom_gas,3))
      allocate(rad_GAS(n_atom_gas))
      if(n_atom_gas==5 )mol_atoms = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      print*,'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'
!      mol_atoms = (/ iSilicon,iOxygenH,iOxygenH /)
      allocate( mol_vec(3,n_atom_GAS) )
! mol_vec(:,1) = 0.0_wp
! mol_vec(:,2) = (/ 0.0_wp, 0.0_wp, ASiO /)
! mol_vec(:,3) = (/ 0.0_wp, -(2.0_wp*((1.0_wp/3.0_wp)*ASiO))*sqrt(2.0_wp), (1.0_wp/3.0_wp)*ASiO /)
      if (n_atom_GAS > 0) call tetra_coords5(ASiO, mol_vec)

print '(3f12.6)',mol_vec
print *,mol_atoms
print '(a)',atom_name(mol_atoms)
!     Initialize the gas molecules
      do i = 1,n_GAS
         !GAS_charge(i,:) =
         GAS_atom(i,:) = mol_atoms(1:n_atom_gas)
         do j = 1, n_atom_gas
            GAS_mass(i,j,:) = amass(mol_atoms(j))
         end do
      end do
      GAS_massi = 1.0_wp/GAS_mass

      rad_GAS = sigLJ_2(mol_atoms)



      allocate(mol_atoms_h2o(n_atom_h2o))
      allocate(h2o_atom(1:n_h2o,n_atom_h2o))
      allocate(h2o_charge(n_h2o,n_atom_h2o))
      allocate(h2o_mass(1:n_h2o,n_atom_h2o,3),h2o_massi(1:n_h2o,n_atom_h2o,3))
      allocate(rad_h2o(n_atom_h2o))
      if (n_atom_h2o == 1) mol_atoms_h2o = (/ iOw /)
      allocate( mol_vec_h2o(3,n_atom_h2o) )
      if (n_atom_h2o == 1) mol_vec_h2o(:,1) = 0.0_wp
! mol_vec_h2o(:,2) = (/ 0.0_wp, 0.0_wp, ??? /)
! mol_vec_h2o(:,3) = (/ 0.0_wp, ???, ??? /)
!if (n_atom_h2o > 0) call h2o_coords(???,???)
print *,'H2O'
print '(3f12.6)',mol_vec_h2o
print *,'mol_atoms_h2o',mol_atoms_h2o
!print '(a)',atom_name(mol_atoms_h2o)
!     Initialize the h2o molecules
      do i = 1,n_h2o
         !h2o_charge(i,:) =
         h2o_atom(i,:) = mol_atoms_h2o(1:n_atom_h2o)
         do j = 1, n_atom_h2o
            h2o_mass(i,j,:) = amass(mol_atoms_h2o(j))
         end do
      end do
if(n_atom_h2o==1) h2o_mass = mOx+2*mH
      h2o_massi = 1.0_wp/h2o_mass

      rad_h2o = sigLJ_2(mol_atoms_h2o)

print*,'abccc'
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
!     nat = 2
!
      allocate(h2o_xyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(r_h2o_uc(n_h2o,n_atom_h2o,3))
      allocate(dr_h2o(n_h2o,n_atom_h2o,3))
      allocate(rcm_h2o_uc(n_h2o,3))
      allocate(vcm_h2o(n_h2o,3))
      allocate(image_h2o_xyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(h2o_fxyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(h2o_vxyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
!
      allocate(GAS_xyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(r_GAS_uc(n_GAS,n_atom_gas,3))
      allocate(dr_GAS(n_GAS,n_atom_gas,3))
      allocate(rcm_GAS_uc(n_GAS,3))
      allocate(vcm_GAS(n_GAS,3))
      allocate(image_GAS_xyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(GAS_fxyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(GAS_vxyz(1:n_GAS,1:n_atom_gas,1:NDIM))
!
      allocate(Ox_xyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(r_O2_uc(n_O2,n_atom_o2,3))
      allocate(dr_O2(n_O2,n_atom_o2,3))
      allocate(rcm_O2_uc(n_O2,3))
      allocate(vcm_O2(n_O2,3))
      allocate(image_Ox_xyz(1:n_O2,1:NDIM,1:NDIM))
      allocate(Ox_fxyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(Ox_vxyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(Ox_atom(1:n_O2,3))
      allocate(Ox_gas_charge(n_O2,3))
      allocate(Ox_mass(1:n_O2,n_atom_o2,3))
      allocate(Ox_massi(1:n_O2,n_atom_o2,3))
!
      allocate(N2_xyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(r_N2_uc(n_N2,n_atom_n2,3))
      allocate(dr_N2(n_N2,n_atom_n2,3))
      allocate(rcm_N2_uc(n_N2,3))
      allocate(vcm_N2(n_N2,3))
      allocate(image_N2_xyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_fxyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_vxyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_atom(1:n_N2,3))
      allocate(N2_gas_charge(n_N2,3))
      allocate(N2_mass(1:n_N2,n_atom_n2,3))
      allocate(N2_massi(1:n_N2,n_atom_n2,3))
!
      allocate(CO2_xyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(r_CO2_uc(n_CO2,n_atom_co2,3))
      allocate(dr_CO2(n_CO2,n_atom_co2,3))
      allocate(rcm_CO2_uc(n_CO2,3))
      allocate(vcm_CO2(n_CO2,3))
      allocate(image_CO2_xyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_fxyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_vxyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_atom(1:n_CO2,n_atom_co2))
      allocate(CO2_gas_charge(n_CO2,n_atom_co2))
      allocate(CO2_mass(1:n_CO2,n_atom_co2,3))
      allocate(CO2_massi(1:n_CO2,n_atom_co2,3))

      allocate( ULJEL(max(n_o2,n_n2,n_co2,n_gas,n_h2o)) )
      allocate( e_bond_CO2(n_co2),e_angle_CO2(n_co2),e_bond_gas(n_gas),e_angle_gas(n_gas) )

! Read in the initial molecular system file header
      select case(ftyp)
      case('xyz')
         call read_xyz_header(io_config,natom,boxl)
         call get_free_file_unit(io_topol)
         open(unit=io_topol,file=trim(coordfile(1:nc - 3))//'con')
      case('pdb')
         call read_pdb_header(io_config,natom,boxl)
 print*,'BOXL==',boxl(1),boxl(2),boxl(3)
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
      

!
      natom_max = natom
!      boxl = boxl/Angstrom
!      boxl2 = boxl/2.0_wp
!      boxli = 1.0_wp/boxl
      call SET_BOX_SIZE(boxl)
!      
      
!rcutel = boxl2
!rcutel2 = boxl2**2
!
      allocate(rxyz(1:natom_max,3),coprxyz(1:natom_max,3))
      allocate(fxyz(1:natom_max,3),vxyz(1:natom_max,3))
      allocate(atom(1:natom_max),latom(1:natom_max),atom2(1:natom_max))
      atom = 0
      allocate(charge(1:natom))
      allocate(proximity(1:natom_max,4))
      proximity = 0
      allocate(mass(1:natom_max,3),massi(1:natom_max,3))


! read in coordinates and atom types
      select case(ftyp)
      case('xyz')
         call read_xyz(io_config,natom,atom,rxyz)
         close(io_config)
      case('pdb')
         call read_pdb(io_config,natom,atom,rxyz)
      end select

      
!      rxyz(1:natom,:) = rxyz(1:natom,:)/Angstrom
!print*,'rxyz(1:natom,:)',rxyz(1:natom,:)

      do i = 1,natom
         mass(i,:) = amass(atom(i))
      end do
      massi(1:natom,:) = 1.0_wp/mass(1:natom,:)

! read in connectivity
      select case(ftyp)
      case('xyz')
         call read_conn(io_topol,proximity)
         close(io_topol)
      case('pdb')
         call read_conect_pdb(io_config,natom,proximity)
         close(io_config)
      end select

      
dLflexible = -5.0_wp/Angstrom
      dLkeep = min(minval(rxyz(:,3)),-15.0_wp/angstrom)

      print *,' dLflexible = ', dLflexible*Angstrom
      print *,'     dLkeep = ', dLkeep*Angstrom

 ! Make a list of the fixed atoms
      call new_list(fixed_list,natom)
      do i = 1,natom
         if (rxyz(i,3) < dLflexible) call add_to_list(fixed_list, i)
      end do
      ! and renumber the fixed atoms to go from 1,2,3,...,nfixed
      nfixed = fixed_list%n

      do i = 1,nfixed
         call swap_atom(i,fixed_list%ind(i))
         !print *,'swap ',i,fixed_list%ind(i)
      end do
      if (nfixed /= count( rxyz(:,3) < dLflexible )) stop 'nfixed is wrong'
      call delete_list(fixed_list) ! finished with this list

      

   print *,'nfixed = ',nfixed
   print *,' natom = ',natom
      boxz = 84.0_wp ! max hight of the layer
      volbox = boxl(1)*boxl(2)*boxl(3)*(Angstrom**3) ! Box volume in A^3
      volbox1 = boxl(1)*boxl(2)*boxz*(Angstrom**3)
      Rbox_min = (/ -boxl2(1), -boxl2(2), dLKeep /)
      Rbox_max = (/  boxl2(1),  boxl2(2), boxl2(3) /)

print '(a,3f12.6)','Rbox_min = ',Rbox_min*Angstrom
print '(a,3f12.6)','Rbox_max = ',Rbox_max*Angstrom
print *,'volbox =',volbox
print *,'volbox1 =',volbox 

      
print*,'abbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'      
      call assign_charge(1,natom)
      call LJ_INIT

      write(*,*) 'sum charge = ',sum(charge(1:natom))
      write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)

      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      n_oc = count(atom(1:natom) == iOxygenC)
      n_ch2 = count(atom(1:natom) == iCarbonH2)
      n_ch3 = count(atom(1:natom) == iCarbonH3)
      ntot = n_si + n_ob + n_oh + n_oc + n_ch3 + n_ch2
      write(*,*) 'ntot = ',n_si + n_ob + n_oh + n_oc + n_ch3 + n_ch2  
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) '# Oc = ',n_oc
      write(*,*) '# Ch2 = ',n_ch2
      write(*,*) '# Ch3 = ',n_ch3
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh+ n_oc + n_ch3 + n_ch2)/ntot
      !write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      fdens_wt = frac_dens_wt()
      write(*,*) 'density % = ',100.0*fdens_wt
      write(*,*) 'N_Si = ',n_si/volbox,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob/volbox,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh/volbox,' Angstrom^-3'
      write(*,*) 'N_Oc = ',n_oc/volbox,' Angstrom^-3'
      write(*,*) 'N_ch2= ',n_ch2/volbox,' Angstrom^-3'
      write(*,*) 'N_ch3 = ',n_ch3/volbox,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh+n_oc*mOx+n_ch2*mch2+n_ch3*mch3)*(1e24_wp/(volbox1*NAvo))
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'



call check_proximity(1,natom,ok,i)
if (.NOT.ok) then
   write(*,*) 'proximity array is NOT consistent'
   write(*,*) 'check_proximity ',ok,i
end if
call check_proximity2(1,natom,atom(1:natom),proximity,ok,i)
if (.NOT.ok) then
  write(*,*) 'proximity array is NOT consistent'
  write(*,*) 'check_proximity2 ',ok,i
  stop
end if


!      do i = nfixed+1,natom
!         if (ALL(proximity(i,:)==0)) then
!         !if (count(proximity(i,:) > 0) == 0) then
!            print *,'free atom: ',i,atom_name(atom(i)),count(proximity(i,:) > 0)
!         end if
!         select case(atom(i))
!         case(iOxygen)
!            if (count(proximity(i,:) > 0) < 2) then
!               print *,i,atom_name(atom(i)),count(proximity(i,:) > 0)
!            end if
!         case(iOxygenH)
!         case(iSilicon)
!            if (count(proximity(i,:) > 0) < ncmax(atom(i))) then
!               print *,i,atom_name(atom(i)),count(proximity(i,:) > 0)
!            end if
!            nc = count( atom(proximity(i,:))==iOxygen )
!            if (nc < ncmax(atom(i))) then
!               print *,i,atom_name(atom(i)),nc,noh(i),nc+noh(i),count(proximity(i,:) > 0)
!               print '(i5,2x,a2,2x,4i6,4(1x,a2))',i,atom_name(atom(i)),proximity(i,:),atom_name(atom(proximity(i,:)))
!            end if
!         case default
!            print *,i,atom(i)
!            stop 'unknown atom type'
!         end select
!      end do
!
!call new_frame_file(imve,'frame',30000)
!call write_box_frame(imve,n_CO2,print_all = .true.)
!call write_frame(imve,1,natom)
!! remove the free oxygens
!         j = natom
!         noxdel = 0
!         do
!            if (atom(j) == iOxygen .or. atom(j) == iOxygenH) then
!            if (count(proximity(j,:) > 0) < 1) then
!               call delete_atom(j)
!               noxdel = noxdel + 1
!            end if
!            end if
!            j = j - 1
!            if (j == 0) exit
!         end do
!! convert singly bonded oxygens to OH
!         do i = 1,natom
!            if (atom(i) == iOxygen) then
!            if (count(proximity(i,:) > 0) < 2) then
!               atom(i) = iOxygenH
!            end if
!            end if
!         end do

!stop

 call new_grid3d(Rbox_min,Rbox_max,del_bin,Grid)

     !print '(a,3f12.6)','Grid cell sides = ',Grid%del
     !print '(a,3i7)',   'Grid cell bins  = ',Grid%nbin
      Nx = Grid%nbin(1)
      Ny = Grid%nbin(2)
      Nz = Grid%nbin(3)

print*,'Nx=',Nx,'Ny=',Ny,'Nz=',Nz
      nbinmax = Nx*Ny*Nz
      volbin = Grid%volbin
print*,'nbinmax =',nbinmax 
      allocate(nden(nx,ny,nz,ntyp),ndenz(Nz,ntyp))
      allocate(HenryK(nx,ny,nz,n_vapor_species))
      allocate(Accessible(nx,ny,nz,n_vapor_species))
      allocate(Lattice(nx,ny,nz,n_vapor_species))

      do i = 1,n_vapor_species
         call new_list(NumAcc(i),Nz)
      end do

!!
      call INIT_NLIST(Rbox_min,Rbox_max,del_nlist,natom_max)
      
      call NEW_NLIST(1,natom)
      call Init_Verlet_List(rv0=5.0_wp/Angstrom,rcut0=2.6_wp/Angstrom,nat_max=natom_max)
      call new_vlist(1,natom)

!      
      call Init_HKNonLattice(natom_max)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
     write(*,*) 'n_cluster = ',n_cluster

!PERCOLATION ANALYSIS
     top_atom = maxval(rxyz(1:natom,3))
print *,'top_atom',top_atom
print *,'Grid%rlower==',Grid%rlower(3),'Grid%deli',Grid%deli(3)
      nz_top_atom = int((top_atom-Grid%rlower(3))*Grid%deli(3)) + 1
print *,'nz_top_atom =',nz_top_atom
      nz_top = nz_top_atom + 3
print *,'nz_top =',nz_top
      z_top = nz_top*Grid%del(3)
print *,'z_top =',z_top
      Rbox_top = (/ boxl2(1), boxl2(2), z_top /)
print *,'Rbox_top =',Rbox_top
print *,'z_top = ',z_top*Angstrom

      call den_profile_grid(nfixed+1,natom,Grid,nden,ndenz)
print*,'denprofgrid'
      Accessible = .true.
      call Henrys_law_grid(probe_CO2,Accessible(:,:,1:Nz_top,iCO2),grid,Nz_top, &
                           n_try_henry,fover_henry,rcutH,henryk(:,:,1:Nz_top,iCO2))
      print*,'henryk(:,:,1:Nz_top,iCO2=',henryk(:,:,1:Nz_top,iCO2)
      where (henryk(:,:,1:nz_top,:) >= void_crit)
             Lattice(:,:,1:nz_top,:) = 1
      elsewhere
            Lattice(:,:,1:nz_top,:) = 0
      end where 


     do i = 1, n_vapor_species 
         ! Use the Hoshen-Kopelman algorithm to label the clusters of connected cells for
         ! each vapor phase reactant
         call hk_lattice_uf( Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i) )
         ! Alternative versions of algorithm (for comparison)
         ! call hk_lattice_3d(Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i))
         ! call hk_lattice_3dnl(Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i))
         !
         ! Next mark the cell which are accessible (from the vapor phase)
         call ACCESSIBLE_REGIONS_Z( Lattice(:,:,1:Nz_top,i), 1, Nz_top, &
                                 Accessible(:,:,1:Nz_top,i), NumAcc(i) )
!print *,i,ncluster(i)
      end do
!print*,'Accessible(:,:,1:Nz_top,iSiOC2H54)=',Accessible(:,:,1:Nz_top,iSiOC2H54)
print*,'Accessible(:,:,1:Nz_top,iCO2)=',Accessible(:,:,1:Nz_top,iCO2)


      do j = 1, n_vapor_species
         print *,'# NumAcc%n = ',NumAcc(j)%n
         do i = 1,NumAcc(j)%n
            print *,i,NumAcc(j)%ind(i)
         end do
      end do

print *,' percolation check  '

! FINISHED PERCOLATION ANALYSIS
!
      do i = 1,n_O2
        Ox_gas_charge(i,1) = q_O_O2
        Ox_gas_charge(i,2) = q_O_O2
        Ox_gas_charge(i,3) = -2.0_wp*q_O_O2
        Ox_atom(i,1) = iO_O2
        Ox_atom(i,2) = iO_O2
        Ox_mass(i,:,:) = mOx
      end do
      Ox_massi = 1.0_wp/Ox_mass
      do i = 1,n_N2
        N2_gas_charge(i,1) = q_n_N2
        N2_gas_charge(i,2) = q_n_N2
        N2_gas_charge(i,3) = -2.0_wp*q_n_N2
        N2_atom(i,1) = in_N2
        N2_atom(i,2) = in_N2
        N2_mass(i,:,:) = mN
      end do
      N2_massi = 1.0_wp/N2_mass
      do i = 1,n_CO2
        CO2_gas_charge(i,1) = q_C_CO2
        CO2_gas_charge(i,2) = q_O_CO2
        CO2_gas_charge(i,3) = q_O_CO2
        CO2_atom(i,1) = iC_CO2
        CO2_atom(i,2) = iO_CO2
        CO2_atom(i,3) = iO_CO2
        CO2_mass(i,1,:) = mC
        CO2_mass(i,2:3,:) = mOx
      end do
      CO2_massi = 1.0_wp/CO2_mass

! Insertion
      !zz = rxyz(:,3)
      zmax = maxval(rxyz(:,3)); print *,'zmax = ',zmax
      zmax1 = (zmax - 4.0_wp) ; print *,'zmax1 = ',zmax1
      zmin = minval(rxyz(:,3)); print *,'zmin = ',zmin
      zmin1 = (zmin + 3.0_wp); print *,'zmin1 = ',zmin1

      Rbox_min1 = (/ -boxl2(1), -boxl2(2),zmin1 /)
      Rbox_max1 = (/  boxl2(1),  boxl2(2),zmax1 /)
!rl = -boxl2
!ru =  boxl2
print*,'Rbox_min1',Rbox_min1
print*,'Rbox_max1',Rbox_max1

!      call insert_mol_GAS(n_h2o,n_atom_h2o,mol_vec_h2o,h2o_xyz,100000,rad_h2o,rl,ru,fover,check_self_overlap) ; print *,'n_h2o = ',n_h2o

!      call insert_mol_GAS(n_GAS,n_atom_gas,mol_vec,GAS_xyz,100000,rad_GAS,zmax,ru,fover,check_self_overlap) ; print *,'n_gas = ',n_gas
      if (allocated(mol_vec)) deallocate(mol_vec)
      if (.not.allocated(mol_vec)) allocate(mol_vec(3,2))

      if (n_O2 > 0) then
      mol_vec(1,1) =  AOO*0.5_wp
      mol_vec(1,2) = -AOO*0.5_wp
!      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_O2,n_atom_o2,mol_vec(:,1:n_atom_o2),Ox_xyz,1000000,rad_Ox,Rbox_min1,Rbox_max1,fover,check_self_overlap)
      print *,'n_o2 = ',n_o2
      end if
      if (n_N2 > 0) then
      mol_vec(1,1) =  ANN*0.5_wp
      mol_vec(1,2) = -ANN*0.5_wp
!      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_N2,n_atom_n2,mol_vec,N2_xyz,1000000,rad_N2,Rbox_min1,Rbox_max1,fover,check_self_overlap)
      print *,'n_n2 = ',n_n2
      end if
      if (allocated(mol_vec)) deallocate(mol_vec)
      if (.not.allocated(mol_vec)) allocate(mol_vec(3,3))
      if (n_CO2 > 0) then
      mol_vec(1,1) = 0.0_wp
      mol_vec(1,2) =  ACO2
      mol_vec(1,3) = -ACO2
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_CO2,n_atom_co2,mol_vec,CO2_xyz,1000000,rad_CO2,Rbox_min1,Rbox_max1,fover,check_self_overlap)
      print *,'n_co2 = ',n_co2
      end if
!
      if (restart) then
         call read_box_frame(io_restart_xyz,nmol = n_GAS,print_all = .true.)
         call read_vel(io_restart_v,nmol = n_GAS,print_all = .true.)
      end if
!
! Set up the uncorrected coords
      do i = 1, n_h2o
         r_h2o_uc(i,1,:) = h2o_xyz(i,1,:)
         do j = 2,n_atom_h2o
            rr = h2o_xyz(i,j,:) - h2o_xyz(i,1,:)
            call pbc(rr)
            r_h2o_uc(i,j,:) =  r_h2o_uc(i,1,:) + rr
         end do
      end do

      do i = 1, n_GAS
         r_GAS_uc(i,1,:) = GAS_xyz(i,1,:)
         do j = 2,n_atom_GAS
            rr = GAS_xyz(i,j,:) - GAS_xyz(i,1,:)
            call pbc(rr)
            r_GAS_uc(i,j,:) =  r_GAS_uc(i,1,:) + rr
         end do
      end do
      do i = 1, n_O2
         r_O2_uc(i,1,:) = Ox_xyz(i,1,:)
         rr = Ox_xyz(i,2,:) - Ox_xyz(i,1,:)
         call pbc(rr)
         r_O2_uc(i,2,:) = r_O2_uc(i,1,:) + rr
      end do
      do i = 1, n_N2
         r_n2_uc(i,1,:) = N2_xyz(i,1,:)
         rr = N2_xyz(i,2,:) - N2_xyz(i,1,:)
         call pbc(rr)
         r_N2_uc(i,2,:) =  r_n2_uc(i,1,:) + rr
      end do
      do i = 1, n_CO2
         r_CO2_uc(i,1,:) = CO2_xyz(i,1,:)
         rr = CO2_xyz(i,2,:) - CO2_xyz(i,1,:)
         call pbc(rr)
         r_CO2_uc(i,2,:) =  r_CO2_uc(i,1,:) + rr
         !
         rr = CO2_xyz(i,3,:) - CO2_xyz(i,1,:)
         call pbc(rr)
         r_CO2_uc(i,3,:) = r_CO2_uc(i,1,:) + rr
      end do


      if (make_movie) then ! write out first frame of movie
         iframe = 0
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_CO2,print_all = .true.)
      end if
!========================================================================
!
      

      call comvel(natom,T_kelvin,vxyz,mass)
      call comvel(n_h2o,n_atom_h2o,T_kelvin,h2o_vxyz,h2o_mass)
      call comvel(n_GAS,n_atom_GAS,T_kelvin,GAS_vxyz,GAS_mass)
      call comvel(n_O2,2,T_kelvin,Ox_vxyz,Ox_mass)
      call comvel(n_N2,2,T_kelvin,N2_vxyz,N2_mass)
      call comvel(n_CO2,3,T_kelvin,CO2_vxyz,CO2_mass)
!
      print '(a/,3g18.8)','net GAS momentum', &
         sum(GAS_vxyz(:,:,1)*GAS_mass(:,:,1)),&
         sum(GAS_vxyz(:,:,2)*GAS_mass(:,:,2)),&
         sum(GAS_vxyz(:,:,3)*GAS_mass(:,:,3))
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      do j = 1,n_atom_GAS
      do i = 1,n_GAS
         SUMX = SUMX + GAS_vxyz(i,j,1)*GAS_mass(i,j,1)
         SUMY = SUMY + GAS_vxyz(i,j,2)*GAS_mass(i,j,2)
         SUMZ = SUMZ + GAS_vxyz(i,j,3)*GAS_mass(i,j,3)
      end do
      end do
      print '(3g18.8)',SUMX,SUMY,SUMZ
      print *

      print '(a/,3g18.8)','net H2O momentum', &
         sum(h2o_vxyz(:,:,1)*h2o_mass(:,:,1)),&
         sum(h2o_vxyz(:,:,2)*h2o_mass(:,:,2)),&
         sum(h2o_vxyz(:,:,3)*h2o_mass(:,:,3))

      print '(a/,3g18.8)','net CO2 momentum', &
         sum(CO2_vxyz(:,:,1)*CO2_mass(:,:,1)),&
         sum(CO2_vxyz(:,:,2)*CO2_mass(:,:,2)),&
         sum(CO2_vxyz(:,:,3)*CO2_mass(:,:,3))
!print *,'CO2_mass(:,1,1:NDIM))',sum(CO2_mass(:,1,1:NDIM))/(n_CO2*3)
!print *,'CO2_mass(:,2,1:NDIM))',sum(CO2_mass(:,2,1:NDIM))/(n_CO2*3)
!print *,'CO2_mass(:,3,1:NDIM))',sum(CO2_mass(:,3,1:NDIM))/(n_CO2*3)

      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      do j = 1,3
      do i = 1,n_CO2
         SUMX = SUMX + CO2_vxyz(i,j,1)*CO2_mass(i,j,1)
         SUMY = SUMY + CO2_vxyz(i,j,2)*CO2_mass(i,j,2)
         SUMZ = SUMZ + CO2_vxyz(i,j,3)*CO2_mass(i,j,3)
      end do
      end do
      print '(3g18.8)',SUMX,SUMY,SUMZ
      print *

      print '(a/,3g18.8)','net solid momentum', &
         sum(vxyz(:,1)*mass(:,1)),&
         sum(vxyz(:,2)*mass(:,2)),&
         sum(vxyz(:,3)*mass(:,3))
      print '(a/,3g18.8)','net solid momentum', &
         dot_product(vxyz(:,1),mass(:,1)),&
         dot_product(vxyz(:,2),mass(:,2)),&
         dot_product(vxyz(:,3),mass(:,3))

!      call voidage_calc(probe_O2,-boxl2,boxl2,100000,voidage(1))
!      call voidage_calc(probe_N2,-boxl2,boxl2,100000,voidage(2))
!      call voidage_calc(probe_CO2,-boxl2,boxl2,100000,voidage(3))
!      call voidage_calc(probe_SiO4,-boxl2,boxl2,100000,voidage(4))
!      call voidage_calc(probe_tip3p,-boxl2,boxl2,100000,voidage(5))
!  !call voidage_profile(probe_O2,-boxl2,boxl2,10000,5.0_wp/Angstrom,nbinv,voidage)
!!print *,'nbin = ',nbinv
!print '(i6,8g18.8)',natom,fdens_wt,dens_gcm3,voidage(1:5)
!open(unit=193, file='voidage.out',status='unknown',position='append')
!write(193,'(i6,8g18.8)')natom,fdens_wt,dens_gcm3,voidage(1:5)
!
!stop



! Rescale each gas to the set temperature
      Ek = kinetic_energy(natom,vxyz,mass)
      Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
      Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
      Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
      Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
      Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

      if (n_h2o > 0) then
         Ndf = 3*(n_atom_h2o*n_h2o)
         T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [h2o] = ',T_kinetic
         h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [h2o] = ',T_kinetic, 'rescaled'
      end if

      if (n_GAS > 0) then
         Ndf = 3*(n_atom_gas*n_GAS)
         T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Gas] = ',T_kinetic
         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Gas] = ',T_kinetic, 'rescaled'
      end if

      if (n_co2 > 0) then
         Ndf = 3*(3*n_CO2)
         T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [CO2] = ',T_kinetic
         CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)
         T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [CO2] = ',T_kinetic, 'rescaled'
      end if

      if (n_o2 > 0) then
         Ndf = 3*(2*n_O2)
         T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [O2] = ',T_kinetic
         Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [O2] = ',T_kinetic, 'rescaled'
      end if

      if (n_n2 > 0) then
         Ndf = 3*(2*n_N2)
         T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [N2] = ',T_kinetic
         N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [N2] = ',T_kinetic, 'rescaled'
      end if

      if (natom > 0) then
      Ndf = 3*(natom)
      T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
      print *,'T_kinetic [Solid] = ',T_kinetic
      vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
      Ek = kinetic_energy(natom,vxyz,mass)
      T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
      print *,'T_kinetic [Solid] = ',T_kinetic, 'rescaled'
      end if

      if (make_movie) then ! write out first frame of velocities
         iframe = 0
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_CO2,print_all = .false.)
      end if

      call set_bond_angle_lists()
      call new_tcf(n_h2o,ntau_msd,nsamp,msd_h2o,'msd_h2o.out')
      call new_tcf(n_h2o,ntau_vcf,nsamp,vcf_h2o,'vcf_h2o.out')
      call new_tcf(n_GAS,ntau_msd,nsamp,msd_gas,'msd_GAS.out')
      call new_tcf(n_GAS,ntau_vcf,nsamp,vcf_gas,'vcf_GAS.out')
      call new_tcf(n_O2,ntau_msd,nsamp,msd_o2,'msd_O2.out')
      call new_tcf(n_O2,ntau_vcf,nsamp,vcf_o2,'vcf_O2.out')
      call new_tcf(n_N2,ntau_msd,nsamp,msd_n2,'msd_N2.out')
      call new_tcf(n_N2,ntau_vcf,nsamp,vcf_n2,'vcf_N2.out')
      call new_tcf(n_CO2,ntau_msd,nsamp,msd_co2,'msd_CO2.out')
      call new_tcf(n_CO2,ntau_vcf,nsamp,vcf_co2,'vcf_CO2.out')

      time = 0.0_wp

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      h2o_fxyz = 0.0_wp
      GAS_fxyz = 0.0_wp
      energy = 0.0_wp
!
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL) ; energy = energy + ULJEL
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_h2o,n_atom_h2o,h2o_atom,h2o_xyz,h2o_fxyz,ULJEL); print*,'ULJ h2o  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_h2o))
do i = 1,n_h2o
 print *,i,ULJEL(i)
end do
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL); print*,'ULJ GAS  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_GAS))
do i = 1,n_gas
 print *,i,ULJEL(i)
end do

      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL); print*,'ULJ O2  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_O2))
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL); print*,'ULJ N2 ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_N2))
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL); print*,'CO2 GAS ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_CO2))
      call force_keating(ukeating)
      call repul_force2(1,natom,repulenergy)
      call GAS_stretching(n_GAS,n_atom_GAS,GAS_xyz,ASiO,kSiO,e_bond_gas)
      energy = energy + sum(e_bond_gas(1:n_GAS))
      call O2_stretching(e_bond_o2) ; energy = energy + e_bond_o2
      call N2_stretching(e_bond_n2) ; energy = energy + e_bond_n2
      call CO2_stretching(e_bond_co2) ; energy = energy + sum(e_bond_co2(1:n_co2))
      call GAS_angle_bending(n_GAS,n_atom_GAS,GAS_xyz,ctheta = cos_OSiO,ktheta = KOSiO,en=e_angle_gas)
      energy = energy + sum(e_angle_gas(1:n_gas))
      call CO2_angle_bending(e_angle_co2)
      energy = energy + sum(e_angle_co2(1:n_co2))
!      if( n_atom_h2o > 1 )then
!         call GAS_stretching(n_h2o,n_atom_h2o,h2o_xyz,A_h2o,k_ohw,e_bond_h2o)
!         energy = energy + sum(e_bond_h2o(1:n_h2o))
!         call GAS_angle_bending(n_h2o,n_atom_h2o,h2o_xyz,ctheta = ctheta_h2o,ktheta = K_h2o,en=e_angle_h2o)
!         energy = energy + sum(e_angle_h2o(1:n_h2o))
!      end if

print *,'e_bond_o2 ',e_bond_o2
print *,'e_bond_n2 ',e_bond_n2
print *,'e_bond_co2 ',sum(e_bond_co2(1:n_co2))
print *,'e_bond_gas ',sum(e_bond_gas)
print *,'e_angle_gas ',sum(e_angle_gas)
do i = 1,n_gas
 print *,i,e_bond_gas(i),e_angle_gas(i)
end do




!======================================================== main loop
      main_loop: do istep = 1, number_moves

         call cpu_time(tp(1))

         coprxyz = rxyz
         call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)
!
         call gas_vverlet_a(dt,h2o_massi,h2o_xyz,h2o_vxyz,h2o_fxyz,dr_h2o)
         call pbc(h2o_xyz)
         r_h2o_uc = r_h2o_uc + dr_h2o

         call gas_vverlet_a(dt,GAS_massi,GAS_xyz,GAS_vxyz,GAS_fxyz,dr_GAS)
         call pbc(GAS_xyz)
         r_GAS_uc = r_GAS_uc + dr_GAS

         call gas_vverlet_a(dt,Ox_massi,Ox_xyz,Ox_vxyz,Ox_fxyz,dr_O2)
         call pbc(Ox_xyz)
         r_O2_uc = r_O2_uc + dr_O2

         call gas_vverlet_a(dt,N2_massi,N2_xyz,N2_vxyz,N2_fxyz,dr_N2)
         call pbc(N2_xyz)
         r_N2_uc = r_N2_uc + dr_N2

         call gas_vverlet_a(dt,CO2_massi,CO2_xyz,CO2_vxyz,CO2_fxyz,dr_CO2)
         call pbc(CO2_xyz)
         r_CO2_uc = r_CO2_uc + dr_CO2

         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:NDIM)) /= CELL(coprxyz(i,1:NDIM))) then
               call NEW_NLIST(1,natom)
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         h2o_fxyz = 0.0_wp
         GAS_fxyz = 0.0_wp
         Ox_fxyz = 0.0_wp
         N2_fxyz = 0.0_wp
         CO2_fxyz = 0.0_wp

         energy = 0.0_wp                        ; call cpu_time(tp(2)); sumtp(1) = sumtp(1) + tp(2)-tp(1)
         call force_keating(ukeating)
         call repul_force2(1,natom,repulenergy)
         call GAS_stretching(n_GAS,n_atom_GAS,GAS_xyz,ASiO,kSiO,e_bond_gas)
         energy = energy + sum(e_bond_gas(1:n_GAS))
         call O2_stretching(e_bond_o2) ; energy = energy + e_bond_o2
         call N2_stretching(e_bond_n2) ; energy = energy + e_bond_n2
         call CO2_stretching(e_bond_co2) ; energy = energy + sum(e_bond_co2(1:n_co2))
         call GAS_angle_bending(n_GAS,n_atom_GAS,GAS_xyz,ctheta = cos_OSiO,ktheta = KOSiO,en=e_angle_gas)
         energy = energy + sum(e_angle_gas(1:n_gas))
         call CO2_angle_bending(e_angle_co2)
         energy = energy + sum(e_angle_co2(1:n_co2))
!        if( n_atom_h2o > 1 )then
!           call GAS_stretching(n_h2o,n_atom_h2o,h2o_xyz,A_h2o,k_ohw,e_bond_h2o)
!           energy = energy + sum(e_bond_h2o(1:n_h2o))
!           call GAS_angle_bending(n_h2o,n_atom_h2o,h2o_xyz,ctheta = ctheta_h2o,ktheta = K_h2o,en=e_angle_h2o)
!           energy = energy + sum(e_angle_h2o(1:n_h2o))
!        end if

!print *,'istep = ',istep
!print *,'e_bond_o2 = ',e_bond_o2
!print *,'e_bond_n2 = ',e_bond_n2
!print *,'e_bond_co2 = ',sum(e_bond_co2)
!print *,'e_angle_co2 = ',sum(e_angle_co2)
!print *,'e_bond_gas = ',sum(e_bond_gas)
!print *,'e_angle_gas = ',sum(e_angle_gas)


!
!     call FORCE_ENERGY_LJ_EL_GAS(1,natom,1,n_GAS,ULJEL) ; energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL); energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL); energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL); energy = energy + ULJEL
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_GAS,3,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_h2o,n_atom_h2o,h2o_atom,h2o_xyz,h2o_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_h2o))
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_GAS))
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_O2))
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_N2))
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_CO2))
! Electrostatic
!      call FORCE_ENERGY_EL(1,natom,1,n_GAS,3,GAS_charge,GAS_xyz,GAS_fxyz,ULJEL) ;  energy = energy + ULJEL
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)  ; energy = energy + ULJEL
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_N2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)  ; energy = energy + ULJEL
!      call FORCE_ENERGY_EL(1,natom,1,n_CO2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)  ;  energy = energy + ULJEL
!




!!call FORCE_ENERGY_LJ_ATTRACTIVE(listOH,n_OH,atom,rxyz,fxyz,ULJAtt)
!call FORCE_ENERGY_LJ_LIST(n_OH,listOH,atom,rxyz,fxyz,ULJAtt)
!!print *,'ULJ OH = ',ULJAtt
!energy = energy + ULJAtt





      call cpu_time(tp(3)); sumtp(2) = sumtp(2) + tp(3) - tp(2)

!
      call vverlet_b(dt,massi,vxyz,fxyz)
      call gas_vverlet_b(dt,h2o_massi,h2o_vxyz,h2o_fxyz)
      call gas_vverlet_b(dt,GAS_massi,GAS_vxyz,GAS_fxyz)
      call gas_vverlet_b(dt,Ox_massi,Ox_vxyz,Ox_fxyz)
      call gas_vverlet_b(dt,N2_massi,N2_vxyz,N2_fxyz)
      call gas_vverlet_b(dt,CO2_massi,CO2_vxyz,CO2_fxyz)

      time = time + dtreal

      if (istep <= nequib) then
         ! print *,sum(fxyz(1:natom,:),dim=1)
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

         Ek_tot = Ek + Ek_O2 + Ek_N2 + Ek_CO2 + Ek_gas + Ek_h2o

         Ndf = 3*(natom + n_atom_gas*n_GAS  + n_atom_h2o*n_h2o + 2*n_O2 + 2*n_N2 + 3*n_CO2) - 3
         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)
!         print*,'Ndf==',Ndf
!         print*,'T_kinetic =',T_kinetic 
!
         if (mod(istep,nsamp) == 0) then
            write(io_energy,'(8g18.8)') time, EK_tot, energy + ukeating + repulenergy*0.5_wp, &
               EK_tot + energy + ukeating + repulenergy*0.5_wp
            write(io_temp,'(8g18.8)') time, T_kinetic
            call flush(17); call flush(18)
         end if
         !if (T_kinetic == 0.0_wp) T_kinetic = 1.0_wp

!         vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
!         h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
!         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
!         Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
!         N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
!         CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)

         call rem_com_momentum(natom,vxyz,mass)
         call rem_com_momentum(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         call rem_com_momentum(n_GAS,n_atom_GAS,GAS_vxyz,GAS_mass)
         call rem_com_momentum(n_O2,2,Ox_vxyz,Ox_mass)
         call rem_com_momentum(n_N2,2,N2_vxyz,N2_mass)
         call rem_com_momentum(n_CO2,3,CO2_vxyz,CO2_mass)

         if (n_h2o > 0) then
            Ndf = 3*(n_atom_h2o*n_h2o)
            T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
            h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_GAS > 0) then
            Ndf = 3*(n_atom_gas*n_GAS)
            T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
            GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_co2 > 0) then
            Ndf = 3*(3*n_CO2)
            T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
            CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_o2 > 0) then
            Ndf = 3*(2*n_O2)
            T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
            Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_n2 > 0) then
            Ndf = 3*(2*n_N2)
            T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
            N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (natom > 0) then
            Ndf = 3*(natom)
            T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
            vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
         end if

      end if

      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then

         if (n_h2o > 0) then
         do i = 1,n_h2o
            mtmp = 0.0_wp
            rcm_h2o_uc(i,:) = 0.0_wp
            vcm_h2o(i,:) = 0.0_wp
            do j = 1,n_atom_h2o
               rcm_h2o_uc(i,:) = rcm_h2o_uc(i,:) + r_h2o_uc(i,j,:)*amass(h2o_atom(i,j))
               vcm_h2o(i,:) = rcm_h2o_uc(i,:) + h2o_vxyz(i,j,:)*amass(h2o_atom(i,j))
               mtmp = mtmp + amass(h2o_atom(i,j))
            end do
            rcm_h2o_uc(i,:) = rcm_h2o_uc(i,:)/mtmp
            vcm_h2o(i,:) = vcm_h2o(i,:)/mtmp
         end do
         end if

         if (n_gas > 0) then
         do i = 1,n_GAS
            mtmp = 0.0_wp
            rcm_GAS_uc(i,:) = 0.0_wp
            vcm_GAS(i,:) = 0.0_wp
            do j = 1,n_atom_gas
               rcm_GAS_uc(i,:) = rcm_GAS_uc(i,:) + r_GAS_uc(i,j,:)*amass(gas_atom(i,j))
               vcm_GAS(i,:) = rcm_GAS_uc(i,:) + GAS_vxyz(i,j,:)*amass(gas_atom(i,j))
               mtmp = mtmp + amass(gas_atom(i,j))
            end do
            rcm_GAS_uc(i,:) = rcm_GAS_uc(i,:)/mtmp
            vcm_GAS(i,:) = vcm_GAS(i,:)/mtmp
         end do
         end if

         do i = 1,n_O2
            rcm_O2_uc(i,:) = ( r_O2_uc(i,1,:) + r_O2_uc(i,2,:) )*0.5_wp
            vcm_O2(i,:) = ( Ox_vxyz(i,1,:) + Ox_vxyz(i,2,:) )*0.5_wp
         end do
         do i = 1,n_N2
            rcm_N2_uc(i,:) = ( r_N2_uc(i,1,:) + r_N2_uc(i,2,:) )*0.5_wp
            vcm_N2(i,:) = ( N2_vxyz(i,1,:) + N2_vxyz(i,2,:) )*0.5_wp
         end do
         do i = 1,n_CO2
            rcm_CO2_uc(i,:) = ( r_CO2_uc(i,1,:)*mC + r_CO2_uc(i,2,:)*mOx + r_CO2_uc(i,3,:)*mOx )/(mOx + mOx + mC)
            vcm_CO2(i,:) = ( CO2_vxyz(i,1,:)*mC + CO2_vxyz(i,2,:)*mOx + CO2_vxyz(i,3,:)*mOx )/(mOx + mOx + mC)
         end do

         if (n_h2o > 0) then
         call msd_samp(msd_h2o,rcm_h2o_uc)
         call vcf_samp(vcf_h2o,vcm_h2o)
         end if
         if (n_gas > 0) then
         call msd_samp(msd_gas,rcm_GAS_uc)
         call vcf_samp(vcf_gas,vcm_GAS)
         end if
 print*,'iiistep=',istep
         call msd_samp(msd_O2,rcm_O2_uc)
         call vcf_samp(vcf_O2,vcm_O2)
         call msd_samp(msd_N2,rcm_N2_uc)
         call vcf_samp(vcf_N2,vcm_N2)
         call msd_samp(msd_CO2,rcm_CO2_uc)
         call vcf_samp(vcf_CO2,vcm_CO2)
      end if
      end if



      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then
         ! Kinetic energy of each species
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

         Ek_tot = Ek + Ek_O2+ Ek_N2 + Ek_CO2 + Ek_gas + Ek_h2o

         Ndf = 3*(natom + n_atom_gas*n_GAS + n_atom_h2o*n_h2o + 2*n_O2 + 2*n_N2 + 3*n_CO2) - 3

         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)

         write(io_energy,'(8g18.8)') time, Ek_tot, energy + ukeating + repulenergy*0.5_wp, &
            Ek_tot + energy + ukeating + repulenergy*0.5_wp
         write(io_temp,'(9g18.8)') time, t_kinetic, &
            2.0_wp*ek/(k_ev*(3*natom)), &
            2.0_wp*Ek_h2o/(k_ev*(3*(n_atom_h2o*n_h2o))), &
            2.0_wp*Ek_gas/(k_ev*(3*(n_atom_gas*n_gas))), &
            2.0_wp*ek_o2/(k_ev*(3*(n_O2*2))), &
            2.0_wp*ek_n2/(k_ev*(3*(n_N2*2))), &
            2.0_wp*ek_co2/(k_ev*(3*(n_CO2*3)))

         call flush(io_energy);call flush(io_temp)
      end if
      end if

      call cpu_time(tp(4)); sumtp(3) = sumtp(3) + tp(4) - tp(3)
      if (mod(istep,nsamp) == 0) then
         write(*,'(i9,8g18.8)') istep,time,tp(2) - tp(1),tp(3) - tp(2),tp(4) - tp(3), tp(4) - tp(1), tp(4) - tp(0)
      end if

      if (make_movie) then
      if (mod(istep,nfreq_movie) == 0) then
         iframe = iframe + 1
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_CO2,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_CO2,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
      end if
      end if

!===================================================== end main loop
      end do main_loop

      if (n_h2o > 0) then
      call  msd_print(msd_h2o,dtreal)
      call  vcf_print(vcf_h2o,dtreal)
      end if
      if (n_gas > 0) then
      call  msd_print(msd_gas,dtreal)
      call  vcf_print(vcf_gas,dtreal)
      end if
 
      call  msd_print(msd_o2,dtreal)
      call  vcf_print(vcf_o2,dtreal)
      call  msd_print(msd_n2,dtreal)
      call  vcf_print(vcf_n2,dtreal)
      call  msd_print(msd_co2,dtreal)
      call  vcf_print(vcf_co2,dtreal)
      print*,'msd print'

      do i = 1,4
         print *,i,sumtp(i)
      end do
      close(io_energy,status='keep')
      close(io_temp,status='keep')

      call get_free_file_unit(io_config_final)
      open(unit=io_config_final,file='final.xyz')
      call write_xyz(io_config_final,1,natom_max,atom,rxyz)
      close(io_config_final,status='KEEP')

CONTAINS

   FUNCTION frac_dens_wt()
      real(wp):: frac_dens_wt
      real(wp):: mass
      integer:: i
      mass = 0.0_wp
      do i = 1,natom
         select case(atom(i))
         case(iSilicon)
            mass = mass + mSi
         case(iOxygen)
            mass = mass + mOx
         case(iOxygenH)
            mass = mass + mOH
         case(iOxygenC)
            mass = mass + mOx     
         case(iCarbonH2)
            mass = mass + mCh2     
         case(iCarbonH3)
            mass = mass + mCh3        
         case default
            stop 'unknown atom type(frac_dens_Wt)'
         end select
      end do
      frac_dens_wt = mass/(8000*mSi + 16000*mOx)
   END FUNCTION frac_dens_wt

END PROGRAM simbox

