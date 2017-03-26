
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

!!>include 'command_line_mod.f90'

MODULE COMMAND_LINE_MOD
   implicit none
CONTAINS

   FUNCTION command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   END FUNCTION

   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(len=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
      CHARACTER(len=1000) :: TMPVAL
      INTEGER,EXTERNAL :: IARGC
      if (NUMBER < 0) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(STATUS)) STATUS = 1
          RETURN
      else if (NUMBER > IARGC()) then
          if (PRESENT(VALUE )) VALUE  = ' '
          if (PRESENT(LENGTH)) LENGTH = 0
          if (PRESENT(STATUS)) STATUS = 2
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
      if (PRESENT(STATUS)) STATUS = 0
      RETURN
   END SUBROUTINE GET_COMMAND_ARGUMENT

END MODULE

!!>include 'command_arg_mod.f90'


MODULE command_arg_mod
   USE precision_mod
   USE command_line_mod
   implicit none
   !
   ! Gets a command argument and assigns it to a variable.
   ! It checks the status of the argument and writes a record to stdout.
   !
   interface get_com_arg
      module procedure get_com_arg_int, get_com_arg_real, &
                       get_com_arg_logical,get_com_arg_string
   end interface
CONTAINS

   SUBROUTINE get_com_arg_int(ia,var,var_name)
      integer,intent(in):: ia
      integer,intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_real(ia,var,var_name)
      integer,intent(in):: ia
      real(wp),intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_logical(ia,var,var_name)
      integer,intent(in):: ia
      logical,intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_string(ia,var,var_name)
      integer,intent(in):: ia
      character(*),intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      var = ''
      var(1:nc) = carg(1:nc)
      write (*,*) var_name,' = ',trim(var)
   END SUBROUTINE

END MODULE command_arg_mod

!!>include 'Keating_parameters_mod_vonAlfthan.f90'

MODULE Keating_parameters_mod
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
   real(wp),parameter:: kkbond(1:3) = (/ KSiO,KSiSi,KSiO /)
   real(wp),parameter:: aabond(1:3) = (/ ASiO,ASiSi,ASiO /)
   real(wp),parameter:: acosang(0:2) = (/ 1.0_wp,1.0_wp/3.0_wp,1.0_wp /)
END MODULE Keating_parameters_mod

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

!!>include 'global_vars_mod.f90'


MODULE global_vars_mod
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, etot
   integer:: natom,nseed,natom_max
   integer:: nattached, nrelax
END MODULE global_vars_mod

!!>include 'HKNonLattice2.f90'

MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster = 0,n_cluster_old = 0
   integer,allocatable:: atomL(:),atomL_old(:)
CONTAINS

   SUBROUTINE Init_HKNonLattice(natom_max)
      integer,intent(in):: natom_max
      integer:: stat
      if (.not.allocated(atomL)) then
         allocate(atomL(natom_max),atomL_old(natom_max),stat = stat)
         if (stat /= 0) then
            print *,'ERROR: Init_HKNonLattice'
            print *,'ALLOCATION FAILED, stat = ',stat
            stop
         end if
      end if
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
      integer:: i,ii,j,nnmax,iNodeNextL,inn,NodeLPmin,k,m
      integer:: RelabL1(NumberOfNodes),RelabL(NumberOfNodes),RelabLB(NumberOfNodes)
      !
      ! STEPS 1,2 & 3 of AL - Futaisi and Tadeusz Patzek algorithm
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

END MODULE HKNonLattice_mod

!!>include 'seaton_mod2.f90'


MODULE seaton_mod
   USE precision_mod
   USE global_vars_mod, only: angstrom, K_ev, qstar, pi
   implicit none
! atom types
   integer,parameter:: iO_Sil = 0
   integer,parameter:: iSi_Sil = 1
   integer,parameter:: iO_OH_sil = 2
   integer,parameter:: iH_OH_sil = 3
   integer,parameter:: iO_H2O = 4
   integer,parameter:: iH_H2O = 5
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
   real(wp),parameter:: bondl_SiO = 1.61_wp/angstrom
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
! LJ parameters
!   real(wp),parameter:: epsi(0:ntyplj) = (/ &
!                        eps_O_sil, &
!                        eps_Si_sil, &
!                        eps_O_OH_sil, &
!                        eps_H_OH_sil, &
!                        eps_O_H2O, &
!                        eps_H_H2O, &
!                        eps_O_CO2, &
!                        eps_C_CO2, &
!                        eps_N_N2, &
!                        eps_O_O2 /)
!   real(wp),parameter:: sigi(0:ntyplj) = (/ &
!                        sig_O_sil, &
!                        sig_Si_sil, &
!                        sig_O_OH_sil, &
!                        sig_H_OH_sil, &
!                        sig_O_H2O, &
!                        sig_H_H2O, &
!                        sig_O_CO2, &
!                        sig_C_CO2, &
!                        sig_N_N2, &
!                        sig_O_O2 /)
!   real(wp),parameter:: sigi2(0:ntyplj) =sigi*0.5_wp
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

!!>include 'constants_mod2.f90'

MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: angstrom,pi
   USE atom_types_mod
   USE keating_parameters_mod
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
!
   real(wp),parameter:: AOH = bondl_OH
   real(wp),parameter:: aSiOSi = pi*(140.0_wp/180.0_wp)
!  Dreiding
   real(wp),parameter:: KOH = 3035.2_wp    ! ((700 kcal/mol)/A2)   eV nm^-2
   real(wp),parameter:: KSiOH = 4.336_wp   ! ((100 kcal/mol)/rad2)
!
! the following are now taken from Seaton
!  real(wp),parameter:: sigma_Obridge,sigma_OH,sigma_O,sigma_Si,sigma_C,sigma_N
!  real(wp),parameter:: bondl_O2,bondl_N2,bondl_SiO,bondl_CO2
!  real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
! united atom SiOH4
   real(wp),parameter:: r_tet = bondl_SiO + sigma_OH*0.5_wp
!
   real(wp),parameter:: sigma(0:ntyp) = (/ sig_O_sil, &
                                          sig_Si_sil,&
                                            sigma_OH,  &
                                        sig_H_OH_sil,&
                                           sig_O_H2O,&
                                           sig_H_H2O /)
   real(wp),parameter:: sigma_2(0:ntyp) = sigma*0.5_wp
   real(wp),parameter:: epsi(0:ntyp) = (/ eps_O_sil, &
                                         eps_Si_sil, &
                                       eps_O_OH_sil, &
                                       eps_H_OH_sil, &
                                          eps_O_H2O, &
                                          eps_H_H2O /)
   real(wp),parameter:: qi(0:ntyp) = (/ q_O_sil, &
                                       q_Si_sil, &
                                     q_O_OH_sil, &
                                     q_H_OH_sil, &
                                        q_O_H2O, &
                                        q_H_H2O /)
END MODULE constants_mod

!!>include 'atom_types_mod.f90'


MODULE atom_types_mod
   USE precision_mod, only: wp
   implicit none
   integer,parameter:: iOxygen = 0
   integer,parameter:: iSilicon = 1
   integer,parameter:: iOxygenH = 2
   integer,parameter:: iHydrogen = 3
   integer,parameter:: iOw = 4
   integer,parameter:: iHw = 5
   integer,parameter:: ntyp = 5
   character(2),parameter:: atom_name(0:ntyp) = (/ 'O ','Si','OH','H ','Ow','Hw' /)
   integer,parameter:: ncmax(0:ntyp) = (/2,4,2,1,1,2/)
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

END MODULE atom_types_mod

!!>include 'coordinates_mod2.f90'



MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:)
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
      read(iu,*) ctmp,ctmp,ctmp,it
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod

!!>include 'connectivity_mod2.f90'


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
         if (i == ncmax(atom(iat))) then
            print *,'>> ',iat,iold,inew
            STOP 'set_proximity: ERROR'
         end if
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
            if (i1 == 0) then
               i1 = proximity(i,2)
               if (i1 == 0) GOTO 100
            end if
            a1 = atom(i1)
            if (a1 /= iSilicon) GOTO 100
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

   PURE FUNCTION OH_groups_OK(iOH1,iOH2)
      logical:: OH_groups_OK
      integer,intent(in):: iOH1,iOH2
      integer:: iSi1,iSi2
      OH_groups_OK = .FALSE.
      iSi1 = proximity(iOH1,1)
      if (iSi1 /= 0) then
         if (atom(iSi1) /= iSilicon) iSi1 = proximity(iOH1,2)
      else
         iSi1 = proximity(iOH1,2)
      end if
      iSi2 = proximity(iOH2,1)
      if (iSi2 /= 0) then
         if (atom(iSi2) /= iSilicon) iSi2 = proximity(iOH2,2)
      else
         iSi2 = proximity(iOH2,2)
      end if
      if (iSi1 == iSi2) RETURN
      if (.NOT.nearest_neighbor2(iSi1,iSi2)) OH_groups_OK = .TRUE.
   END FUNCTION

END MODULE connectivity_mod

!!>include 'charges_implicitH_mod.f90'

MODULE charges_mod
   USE precision_mod
   USE atom_types_mod
   USE constants_mod
   implicit none
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
   END SUBROUTINE

END MODULE charges_mod

!!>include 'probe_mol_mod2.f90'

MODULE probe_mol_mod
   USE precision_mod
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
      integer,allocatable:: nc(:),con(:,:)
   END TYPE probe_mol
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

   SUBROUTINE WRITE_PROBE_MOL(iu,P)
      USE atom_types_mod
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      integer:: i
      write(iu,*) P%n
      write(iu,*) 'Molecule'
      do i = 1,P%n
         write(iu,'(a2,3(1x,f14.8))') atom_name(P%atom(i)),P%r(1:3,i)
      end do
   END SUBROUTINE

END MODULE probe_mol_mod

!!>include 'probe_mol_init_mod.f90'


MODULE probe_mol_init_mod
   USE precision_mod
   USE constants_mod
   USE probe_mol_mod
   implicit none
!
   type(probe_mol):: probe_SiOH4, probe_SiO4, probe_tip3p,probe_tet, probe_SiO3
   type(probe_mol):: probe_H2O,probe_N2,probe_O2,probe_CO2,probe_tet1
!
CONTAINS

   SUBROUTINE init_probe_mols()
      USE tetra_coords_mod
      USE atom_types_mod
      USE constants_mod
      USE charges_mod
! H2O, United atom
      call new_probe_mol(probe_H2O,1,sigma_H2O*0.5_wp)
! Si(OH)4 United atom
      call new_probe_mol(probe_tet1,1,r_tet)
! H2O
      call new_probe_mol(probe_TIP3P,3)
      call h2o_coords(bondl_H2O,angle_H2O,probe_TIP3P%r)
      probe_TIP3P%atom = (/ iO_H2O, iH_H2O, iH_H2O /)
      probe_TIP3P%rad = sigma_2(probe_TIP3P%atom)
      probe_TIP3P%q = qi(probe_TIP3P%atom)
! Si(OH)4
      call new_probe_mol(probe_SiOH4,9)
      call SiOH4_coords(bondl_SiO,bondl_OH,angle_SiOH,probe_SiOH4%r)
      probe_SiOH4%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH, &
                            iHydrogen,iHydrogen,iHydrogen,iHydrogen /)
      probe_SiOH4%rad = sigma_2(probe_SiOH4%atom)
      probe_SiOH4%q = qi(probe_SiOH4%atom)
      probe_SiOH4%q(1) = qsi(4)
! SiO4, just the 4 OH Oxygens
      call new_probe_mol(probe_SiO4,4,sigma_OH*0.5_wp)
      call tetra_coords(bondl_SiO,probe_SiO4%r)
      probe_SiO4%atom = (/ iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO4%rad = sigma_2(probe_SiO4%atom)
!
! SiO4, Si & 4 OH Oxygens
      call new_probe_mol(probe_tet,5,sigma_OH*0.5_wp)
      call tetra_coords5(bondl_SiO,probe_tet%r)
      probe_tet%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      probe_tet%rad = sigma_2(probe_tet%atom)
      probe_tet%q = qi(probe_tet%atom)
      probe_tet%q(1) = qsi(4)
      allocate(probe_tet%nc(5),probe_tet%con(5,4))
      probe_tet%nc = (/ 4,1,1,1,1/)
      probe_tet%con(1,:) = (/ 2,3,4,5 /)
      probe_tet%con(2,:) = (/ 1,0,0,0 /)
      probe_tet%con(3,:) = (/ 1,0,0,0 /)
      probe_tet%con(4,:) = (/ 1,0,0,0 /)
      probe_tet%con(5,:) = (/ 1,0,0,0 /)
!
! SiO3, Si & 3 OH Oxygens
      call new_probe_mol(probe_SiO3,4,sigma_OH*0.5_wp)
      probe_SiO3%atom = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH /)
      probe_SiO3%rad = sigma_2(probe_SiO3%atom)
      allocate(probe_SiO3%nc(4),probe_SiO3%con(4,4))
      probe_SiO3%nc = (/ 3,1,1,1 /)
      probe_SiO3%con(1,:) = (/ 2,3,4,0 /)
      probe_SiO3%con(2,:) = (/ 1,0,0,0 /)
      probe_SiO3%con(3,:) = (/ 1,0,0,0 /)
      probe_SiO3%con(4,:) = (/ 1,0,0,0 /)
      probe_SiO3%r(:,1) = probe_tet%r(:,1)
      probe_SiO3%r(:,2) = probe_tet%r(:,3)
      probe_SiO3%r(:,3) = probe_tet%r(:,4)
      probe_SiO3%r(:,4) = probe_tet%r(:,5)

! O2
      call new_probe_mol(probe_O2,2,sig_O_O2*0.5_wp)
      probe_O2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_O2*0.5_wp /)
      probe_O2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_O2*0.5_wp /)
! N2
      call new_probe_mol(probe_N2,2,sig_N_N2*0.5_wp)
      probe_N2%r(:,1) = (/ 0.0_wp,0.0_wp, bondl_N2*0.5_wp /)
      probe_N2%r(:,2) = (/ 0.0_wp,0.0_wp, -bondl_N2*0.5_wp /)
! CO2
      call new_probe_mol(probe_CO2,3,sig_O_CO2*0.5_wp)
      probe_CO2%rad(1) = sig_C_CO2*0.5_wp
      probe_CO2%r(:,1) = (/ 0.0_wp,0.0_wp,0.0_wp /)
      probe_CO2%r(:,2) = (/ 0.0_wp,0.0_wp, bondl_CO2 /)
      probe_CO2%r(:,3) = (/ 0.0_wp,0.0_wp, -bondl_CO2 /)
! Si(OH)4
      write(*,*) 'SiOH4 just OH'
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
      write(*,*) 'SiO4 ... just OH'
      write(*,*) probe_SiO4%n
      write(*,'(9f12.6)') probe_SiO4%rad
      write(*,'(3f12.6)') probe_SiO4%r
! Si(OH)4
      write(*,*) 'SiO4'
      write(*,*) probe_tet%n
      write(*,'(9i6)') probe_tet%atom
      write(*,'(9f12.6)') probe_tet%rad
      write(*,'(9f12.6)') probe_tet%q
      write(*,'(3f12.6)') probe_tet%r
! Si(OH)3
      write(*,*) 'SiO3 ... just OH'
      write(*,*) probe_SiO3%n
      write(*,'(9i6)') probe_SiO3%atom
      write(*,'(9f12.6)') probe_SiO3%rad
      write(*,'(3f12.6)') probe_SiO3%r

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

END MODULE probe_mol_init_mod

!!>include 'atom_list_mod.f90'


MODULE atom_list_mod
!
   integer,private,parameter:: n_atoml_max = 10000
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
         if (lst%i(i) == iat) exit
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

!!>include 'bond_angle_types_mod.f90'

MODULE bond_angle_types_mod
    USE precision_mod
    USE atom_types_mod
    USE constants_mod
    implicit none

CONTAINS

   SUBROUTINE bond_type(ia1,ia2,kbond,abond)
      integer,intent(in):: ia1,ia2
      real(wp),intent(out):: kbond,abond
      if (compare_bond(atom(ia1),atom(ia2),iSilicon,iOxygen)) then
         kbond = kSiO
         abond = aSiO
      else if (compare_bond(atom(ia1),atom(ia2),iSilicon,iOxygenH)) then
         kbond = kSiO
         abond = aSiO
      else if (compare_bond(atom(ia1),atom(ia2),iSilicon,iSilicon)) then
         kbond = KSiSi
         abond = ASiSi
         stop 'Si - Si bond not allowed'
      else if (compare_bond(atom(ia1),atom(ia2),iOxygenH,iHydrogen)) then
         kbond = KOH
         abond = AOH
      else
         stop 'unknown bond'
      end if
   END SUBROUTINE bond_type


   SUBROUTINE angle_type(ia1,ia2,ia3,kang,ctheta)
      integer,intent(in):: ia1,ia2,ia3
      real(wp),intent(out):: kang,ctheta
      if ( compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygen, iSilicon,iOxygen) .or.&
           compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygenH,iSilicon,iOxygenH).or.&
           compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygen, iSilicon,iOxygenH) ) then
         kang = KOSiO
         ctheta = -1.0_wp/3.0_wp
      else if (compare_angle(atom(ia1),atom(ia2),atom(ia3),iSilicon,iOxygenH,iHydrogen)) then
         kang = KSiOH
         ctheta = cos(angle_SiOH)
      else if (compare_angle(atom(ia1),atom(ia2),atom(ia3),iSilicon,iOxygen,iSilicon)) then
         kang = KSiOSi
         ctheta = -1.0_wp
      else
         stop 'unknown angle'
      end if
   END SUBROUTINE angle_type


   logical pure FUNCTION compare_bond(i1,i2,typ1,typ2)
      integer,intent(in):: i1,i2,typ1,typ2
      compare_bond = (i1 == typ1.and.i2 == typ2).or. &
                     (i1 == typ2.and.i2 == typ1)
   END FUNCTION compare_bond


   logical pure FUNCTION compare_angle(i1,i2,i3,typ1,typ2,typ3)
      integer,intent(in):: i1,i2,i3,typ1,typ2,typ3
      compare_angle = (i1 == typ1.and.i2 == typ2.and.i3 == typ3).or. &
                      (i1 == typ3.and.i2 == typ2.and.i3 == typ1)
   END FUNCTION compare_angle

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
CONTAINS

   SUBROUTINE bond_list()
      integer:: ia1,ia2,j
      nbondtot = 0
      if (.not.allocated(nbond)) then
         allocate( ibond(2,natom_max*4),nbond(natom_max) )
         allocate( iang(3,natom_max*6) )
         allocate( kbond(natom_max*4),abond(natom_max*4))
         allocate( kang(natom_max*6),ctheta(natom_max*6))
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
               call bond_type(ia1,ia2,kbond(nbondtot),abond(nbondtot))
            end if
         end do
      end do
   END SUBROUTINE bond_list


   SUBROUTINE bond_angle_list()
      integer:: ia1,ia2,ia3,ib,ic
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
            call bond_type(ia1,ia2,kbond(nbondtot),abond(nbondtot))
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
               call angle_type(ia1,ia2,ia3,kang(nang),ctheta(nang))
            end if
         end do
      end do
      end do
   END SUBROUTINE bond_angle_list


   SUBROUTINE atom_bond_angle(ia,nbnd,ibnd,kbnd,abnd,nan,ian,kan,can)
      integer,intent(in):: ia
      integer,intent(out):: nbnd,ibnd(:,:),nan,ian(:,:)
      real(wp),intent(out):: kbnd(:),abnd(:),kan(:),can(:)
      integer:: L,M,j,k
! the bonds involving atom ia
      nbnd = 0
      do L = 1, ncmax(atom(ia))
         j = proximity(ia,L)
         if (j == 0) cycle
         nbnd = nbnd + 1
         ibnd(1,nbnd) = ia
         ibnd(2,nbnd) = j
         call bond_type(ia,j,kbnd(nbnd),abnd(nbnd))
      end do

! the angles with ia as the centre atom
      nan = 0
      do L = 1, ncmax(atom(ia)) - 1
         j = proximity(ia,L)
         if (j == 0) cycle
         do M = L + 1, ncmax(atom(ia))
            k = proximity(ia,M)
            if (k == 0) cycle
            nan = nan + 1
            ian(1,nan) = j
            ian(2,nan) = ia
            ian(3,nan) = k
            call angle_type(j,ia,k,kan(nan),can(nan))
         end do
      end do

! the angles with ia as an end atom
      do L = 1, ncmax(atom(ia))
         j = proximity(ia,L)
         if (j == 0) cycle
         do M = 1, ncmax(atom(j))
            k = proximity(j,M)
            if (k == ia .or. k == 0) cycle
            nan = nan + 1
            ian(1,nan) = ia
            ian(2,nan) = j
            ian(3,nan) = k
            call angle_type(ia,j,k,kan(nan),can(nan))
         end do
      end do
   END SUBROUTINE atom_bond_angle

END MODULE bond_angle_list_mod

!!>include 'frames_mod2.f90'

MODULE FRAMES_MOD
   USE precision_mod, only: wp
   implicit none
   public:: new_frame_file, write_frame

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
      USE bond_angle_list_mod
      USE global_vars_mod, only: angstrom
      USE connectivity_mod
      USE atom_types_mod, only: atom, atom_name, ncmax
      USE coordinates_mod, only: rxyz, boxl2, boxl
      integer,intent(in):: iu,ifirst,ilast
      integer:: i,j,ib,ia,iu2
      logical:: lopen,fix_bonds
      character(len=132):: frame_file
      write(iu,'(a6,2x,i7)') 'REMARK',ilast - ifirst + 1
      write(iu,'(a6,1x,3(1x,f8.4))') 'REMARK',boxl*Angstrom
      do i = ifirst,ilast
         write(iu,110)'HETATM',i,atom_name(atom(i)),'   ',i,rxyz(i,1:3)*angstrom
      end do
      call bond_list
!
      inquire(unit=iu,name=frame_file)
      j = len_trim(frame_file)
      frame_file = frame_file(1:j - 3)//'scr'
      fix_bonds = .false.
      do i = 1,nbondtot
         ia = ibond(1,i)
         ib = ibond(2,i)
         if (ANY( ABS(rxyz(ia,1:3) - rxyz(ib,1:3)) > 10.0_wp/angstrom )) then
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
      write(iu,110)'HETATM',99991,'AR','   ',99991, -boxl2*angstrom, -boxl2*angstrom, -boxl2*angstrom
      write(iu,110)'HETATM',99992,'AR','   ',99992, -boxl2*angstrom, boxl2*angstrom, -boxl2*angstrom
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2*angstrom, -boxl2*angstrom, -boxl2*angstrom
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2*angstrom, boxl2*angstrom, -boxl2*angstrom
      write(iu,110)'HETATM',99995,'AR','   ',99995, -boxl2*angstrom, -boxl2*angstrom,boxl2*angstrom
      write(iu,110)'HETATM',99996,'AR','   ',99996, -boxl2*angstrom, boxl2*angstrom,boxl2*angstrom
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2*angstrom, -boxl2*angstrom,boxl2*angstrom
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2*angstrom, boxl2*angstrom,boxl2*angstrom
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


!!>include 'force_keating.f90'

MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none
CONTAINS

   SUBROUTINE energy_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: rl1(3),r1sq,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         rl1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         ebond = ebond + kbond(ib)*(r1sq - abond(ib))**2
      end do
      ebond = 0.5_wp*ebond
!
!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         eangle = eangle + kang(ia)*(cos_theta - ctheta(ia))**2
      end do
      eangle = 0.5_wp*eangle
   END SUBROUTINE energy_keating


   SUBROUTINE force_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: k,dist,rl1(3),r1sq,ktheta,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         k = KSiO
         dist = ASiO
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - k*(r1sq - dist)*rl1(1:3)/r1sq
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + k*(r1sq - dist)*rl1(1:3)/r1sq
         ebond = ebond + 0.5_wp*k*(r1sq - dist)**2
      end do

!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         coef = ktheta*(cos_theta - ctheta(ia))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
         fxyz(ia2,1:3) = fxyz(ia2,1:3) - coef*(-(rl1 + rl2)/coef2 &
                                            + rl1l2*rl1/coef1 &
                                            + rl1l2*rl2/coef3)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - coef*(rl2/coef2 - rl1l2*rl1/coef1)
         fxyz(ia3,1:3) = fxyz(ia3,1:3) - coef*(rl1/coef2 - rl1l2*rl2/coef3)
         eangle = eangle + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE force_keating

END MODULE force_keating_mod


!!>include 'nlist_mod2.f90'


MODULE NLIST_MOD
      USE precision_mod, only: wp
      USE coordinates_mod, only: rxyz,boxl2
      implicit none
      PRIVATE
      integer,public:: neigh
      real(wp),private:: delx,dely,delz,delmin
      real(wp),private:: delxi,delyi,delzi
      integer:: ncelx,ncely,ncelz,ncelt
      integer,parameter:: ncelmax = 150000,neighmx = 350
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

   SUBROUTINE INIT_NLIST(Lx,Ly,Lz,Rc,natom_max)
      real(wp),intent(in):: Lx,Ly,Lz,Rc
      integer,intent(in):: natom_max
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
      LL = 0
      HOC = 0
      LR = 0
   END SUBROUTINE

   PURE FUNCTION CELL(XYZ)
!     determines cell number for position x,y,z
      real(wp),intent(in)::  XYZ(:)
      integer:: CELL
      CELL = INT((XYZ(1) + boxl2)*delxi) &
           + INT((XYZ(2) + boxl2)*delyi)*ncelx &
           + INT((XYZ(3) + boxl2)*delzi)*ncelx*ncely + 1
      RETURN
   END FUNCTION CELL

   SUBROUTINE NEW_NLIST(ifirst,ilast)
!     makes a new neighbour list using the linked-list algorithm
      integer,intent(in):: ifirst,ilast
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head - of - chain
      ! make linked list:
      do i = ifirst,ilast
         ! determine cell number
         ic = CELL(rxyz(i,:))
! debugging
if (ic < 1 .or.ic > ncelt ) then
   print *,i,rxyz(i,:)
   stop 'NEW_NLIST: ic out of bounds'
end if
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

!!>include 'repul_force_NL_vonAlfthan2.f90'


MODULE repul_force_mod
   USE precision_mod, only: wp
   implicit none
   real(wp):: repulenergy

CONTAINS
   SUBROUTINE repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2,r(3)
      integer:: M,L,ic,nc,inn,nn,neigbors(27)
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repulenergy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      !if (atom(L) == iHydrogen) CYCLE outer_atom_loop
      r = rxyz(L,:)
      ic = CELL(r) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            !if (atom(M) == iHydrogen) GOTO 100
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            rr(1:3) = rxyz(M,1:3) - r(1:3)
            call pbc(rr)
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            if (r2 <= rden2) then
               if (nearest_neighbor2(L,M)) GOTO 100
               repulenergy = repulenergy + (r2 - rden2)**2
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      repulenergy = repulenergy*0.5_wp*ff
   END SUBROUTINE repul_energy

   SUBROUTINE repul_energy2(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom,iHydrogen
      USE connectivity_mod
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp):: r2,dx,dy,dz,x,y,z
      integer:: M,L,ic,nc,inn,nn,neigbors(27),i,j,k
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repulenergy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      !if (atom(L) == iHydrogen) CYCLE outer_atom_loop
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
            !if (atom(M) == iHydrogen) GOTO 100
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            dx = rxyz(M,1) - x
            dy = rxyz(M,2) - y
            dz = rxyz(M,3) - z
            if (abs(dx) > boxl2) dx = dx - sign(boxl,dx)
            if (abs(dy) > boxl2) dy = dy - sign(boxl,dy)
            if (abs(dz) > boxl2) dz = dz - sign(boxl,dz)
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
               repulenergy = repulenergy + 0.5_wp*ff*(r2 - rden2)**2
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
   END SUBROUTINE repul_energy2


   SUBROUTINE repul_force(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(3),r2,r(3)
      integer:: M,L,ic,nc,inn,nn,neigbors(27)
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repulenergy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      !if (atom(L) == iHydrogen) CYCLE outer_atom_loop
      r = rxyz(L,:)
      ic = CELL(r) ! link list
      CALL NEIGCELL(ic,1,nn,neigbors)
      cell_loop: do inn = 1,nn
         nc = neigbors(inn)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            !if (atom(M) == iHydrogen) GOTO 100
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            rr(1:3) = rxyz(M,1:3) - r(1:3)
            call pbc(rr)
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            if (r2 <= rden2) then
               if (nearest_neighbor2(L,M)) GOTO 100
               repulenergy = repulenergy + (r2 - rden2)**2
               fxyz(L,:) = fxyz(L,:) + 2.0_wp*ff*(r2 - rden2)*(rr)
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      repulenergy = repulenergy*0.5_wp*ff
   END SUBROUTINE repul_force

   SUBROUTINE repul_force2(ifirst,ilast)
      USE coordinates_mod
      USE atom_types_mod, only: atom,iHydrogen
      USE connectivity_mod
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp):: r2,dx,dy,dz,x,y,z
      integer:: M,L,ic,nc,inn,nn,neigbors(27),i,j,k
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repulenergy = 0.0_wp
      outer_atom_loop: do L = ifirst,ilast
      !if (atom(L) == iHydrogen) CYCLE outer_atom_loop
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
            !if (atom(M) == iHydrogen) GOTO 100
            !if (nearest_neighbor2(L,M)) GOTO 100
            if (M == L) GOTO 100
            dx = rxyz(M,1) - x
            dy = rxyz(M,2) - y
            dz = rxyz(M,3) - z
            if (abs(dx) > boxl2) dx = dx - sign(boxl,dx)
            if (abs(dy) > boxl2) dy = dy - sign(boxl,dy)
            if (abs(dz) > boxl2) dz = dz - sign(boxl,dz)
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
               repulenergy = repulenergy + 0.5_wp*ff*(r2 - rden2)**2
               fxyz(L,1) = fxyz(L,1) + 2.0_wp*ff*(r2 - rden2)*dx
               fxyz(L,2) = fxyz(L,2) + 2.0_wp*ff*(r2 - rden2)*dy
               fxyz(L,3) = fxyz(L,3) + 2.0_wp*ff*(r2 - rden2)*dz
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
   END SUBROUTINE repul_force2

END MODULE repul_force_mod

!include 'quat2mat_mod.f90'
!include 'ran_point_sphere_mod.f90'
!include 'insert_mod.f90'
!include 'lj_el_mod_O2_N2_CO2.f90'
!include 'force_en_el_mod.f90'



MODULE condensation_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE condensation(iO1,iO2)
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      integer,intent(inout):: iO1,iO2
      integer:: i,ii,ilst(3),itmp
      integer:: iSi2,iSi1
!
!      H          H        O
!      |    +     |  -->  / \  + H2O
!   Si-O       Si-O
!
!
!      |        |           iO1           iO2
!     iO1   +  iO2   -->    /  \   +      /  \
!     /        /           /    \
!   iSi1     iSi2        iSi1  iSi2
!
!! iO2 cannot be a seed
!      if (iO2 <= nseed) then
!         if (iO1 <= nseed) RETURN  ! don't allow 2 seed OH groups to react
!         itmp = iO2
!         iO2 = iO1
!         iO1 = itmp
!      end if

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
!
      atom(iO1) = iOxygen
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO1, 0, iSi2)
      call delete_atom(iO2)
      RETURN
   END SUBROUTINE

END MODULE condensation_mod


MODULE bond_switch_O_mod
   USE precision_mod, only: wp
   IMPLICIT NONE
CONTAINS

   SUBROUTINE bond_switch_O(iOb,success)
      USE atom_types_mod
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      implicit none
      integer,intent(in):: iOb
      logical,intent(out):: success
      integer:: iO1,iO2
      integer:: iSi1,iSi2,iSi3,iSi4
!     logical:: consistent

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

      iSi1 = proximity(iOb,1)
      iSi2 = proximity(iOb,2)
      if (iSi1 == 0) stop 'oxygen malfunction 1'
      if (iSi2 == 0) stop 'oxygen malfunction 2'

! randomly select an oxygen iO1 attached to iSi1 (not iOb)
      do
         iO1 = proximity(iSi1,int(rand()*ncmax(atom(iSi1))) + 1)
         if (iO1 /= iOb) exit
      end do
! repeat for iSi2
      do
         iO2 = proximity(iSi2,int(rand()*ncmax(atom(iSi2))) + 1)
         if (iO2 /= iOb) exit
      end do

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



PROGRAM SIMBOX
      USE precision_mod
      USE command_arg_mod
      USE atom_types_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE nlist_mod
      USE bond_angle_list_mod
      USE rand_mod
      USE HKNonLattice_mod
      USE frames_mod
      USE condensation_mod
      USE bond_switch_O_mod
      USE probe_mol_mod
      USE probe_mol_init_mod
      implicit none
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: NAvo = 6.02214e23_wp
      real(wp):: n_si,n_ob,n_oh
      real(wp):: r3(3),rn1,rn2,rn3
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      integer:: natom_old,natom_old1
      integer,allocatable:: atom_old(:),proximity_old(:,:),proximity2(:,:)
      integer,allocatable:: atom_old1(:),proximity_old1(:,:)
      real(wp),allocatable:: rxyz_old(:,:),bondl(:),rxyz_old1(:,:)
      integer,allocatable:: atomc(:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l,frem,dens_gcm3,volbox,fdens_wt,bondlmax
      real(wp):: r1(3),r2(3),rmax,rmax2,rl1,rl2,rmaxSiSi,rmaxSiSi2
      integer:: narg,len,status,noxdel,iSi1,iSi2
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz,isil,jj,ibmax(1)
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl,lst4(4),nc
      integer:: ipairs(4),ip(3,4),itry(3),ntry,ia,i1,i2,a1,a2
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz,irand_seed0
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster,success
      character(len=132) infile,ctmp,outfile,c6*6,c5*5,ftyp*3
      type(atom_list):: lstSi
!
      ip(1,:) = (/ 1,2,3,4 /)
      ip(2,:) = (/ 1,4,3,2 /)
      ip(3,:) = (/ 1,3,2,4 /)
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_com_arg(0, ctmp, "command")
      if (narg /= 6) then
         write(*,*)'usage :'
         write(*,*) trim(ctmp),' frem rmax rmaxSiSi irand_seed  infile  outfile'
         write(*,*)
         write(*,*) 'frem = fractional density (weight) required'
         write(*,*) "rmax = maximum radius for condensing OH's [Angstrom]"
         write(*,*) 'max radius for Silicas [Angstrom]'
         write(*,*) 'irand_seed = seed for RNG'
         write(*,*) 'infile = name of input file'
         write(*,*) 'outfile = base name of output file'
         stop 'error: incorrect no. of command line arguments'
      end if
      call get_com_arg(1, frem, "fractional density (weight) required")
      call get_com_arg(2, rmax, "max radius for condensing OH's")
      call get_com_arg(3, rmaxSiSi, "max radius for Silicas")
      call get_com_arg(4, irand_seed0, "seed for RNG")
      call get_com_arg(5, infile, "input file")
      call get_com_arg(6, outfile, "base name of output file")
      irand_seed = irand_seed0
      write(*,'(a,2f9.6)') 'Check: first 2 random numbers are ',rand(),rand()
      rmax = rmax/Angstrom
      rmax2 = rmax**2
      rmaxSiSi = rmaxSiSi/Angstrom
      rmaxSiSi2 = rmaxSiSi**2

!
      open(unit=14,file=trim(infile))
      nseed = 0
!
      nc = len_trim(infile)
print *,"'",trim(infile),"'"
print *,nc
      ftyp = infile(nc - 2:nc)
      select case(ftyp)
      case('xyz')
         read (14,*) natom
         read (14,*) boxl
         open(unit=15,file=trim(infile(1:nc - 3))//'con')
      case('pdb')
         read (14,*) c6,natom  ;print *,'natom = ',natom
         read (14,*) c6,boxl   ;print *,' boxl = ',boxl
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
!
      natom_max = natom
      boxl = boxl/Angstrom
      boxl2 = boxl*0.5_wp
      boxli = 1.0_wp/boxl
      volbox = (boxl*boxl*boxl)*(Angstrom**3)
!
      allocate(rxyz(1:natom_max,3),rxyz_old(1:natom_max,3))
      allocate(atom(1:natom_max),atom_old(1:natom_max))
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(cluster_list_x(natom_max),cluster_list_y(natom_max),cluster_list_z(natom_max))
      allocate(proximity(natom_max,4),proximity_old(natom_max,4),proximity2(natom_max,4))
      allocate(rxyz_old1(1:natom_max,3))
      allocate(atom_old1(1:natom_max))
      allocate(proximity_old1(natom_max,4))
!
      select case(ftyp)
      case('xyz')
         do i = 1,natom
            read (14,*) ctmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
         close(14)
      case('pdb')
         do i = 1,natom
            read (14,'(a6)',advance = 'no') c6
            read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
      end select
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom

      proximity = 0
      select case(ftyp)
      case('xyz')
         read(15,*) itmp
         do i = 1,natom
            read(15,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(15)
      case('pdb')
         do i = 1,natom
            read(14,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(14)
      end select
!
! randomly translate the initial structure
!     rn1 = rand()
!     rn2 = rand()
!     rn3 = rand()
!     rxyz(1:natom,1) = rxyz(1:natom,1) + (2.0_wp*rn1 - 1.0_wp)*boxl
!     rxyz(1:natom,2) = rxyz(1:natom,2) + (2.0_wp*rn2 - 1.0_wp)*boxl
!     rxyz(1:natom,3) = rxyz(1:natom,3) + (2.0_wp*rn3 - 1.0_wp)*boxl
!     call pbc(rxyz)
!
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      write(*,*) 'density % = ',100.0*frac_dens_wt()
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
!
! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
         stop 'initially there should be one cluster'
      end if

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
      call check_proximity2(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent (2) ',ii
         stop
      end if

!do i = 1,natom
!if (atom(i) ==iSilicon) then
!   do j =1,4
!      ii = proximity(i,j)
!      if (ii == 0) STOP '0 error in proximity'
!   end do
!end if
!end do
!do i = 1,natom
!if ( all(proximity(i,:) == 0) ) then
!   write(*,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
!   stop 'error in proximity'
!end if
!end do


      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom,natom_max)
      call NEW_NLIST(1,natom)




! Make array of silicon atoms
      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == isilicon) call add_to_list(lstSi,i)
      end do
      print *,'no. of silicons = ',lstSi%n
      print *,count(atom(1:natom) == isilicon)
!
! Perfom the etching procedure
!
      isil = 0
      etching: do
         isil = isil + 1

         ! Make array of silicon atoms
         lstSi%n = 0
         do i = 1,natom
            if (atom(i) == isilicon) call add_to_list(lstSi,i)
         end do

         ! pick a Silicon to remove
100      CONTINUE
         ! print *,'lstSi%n ',lstSi%n
         if (lstSi%n == 0) STOP 'NO suitable Silicon found for removal'
         call rand_from_list(lstSi,iat)
!do
!   iat = int(rand()*natom) + 1
!   if (atom(iat) == iSilicon) exit
!end do
         lst4 = proximity(iat,:)

! store old config
         call store_config()

         where (lst4 == natom) lst4 = iat
         call delete_atom(iat)

         itry = (/ 1,2,3 /)
         ntry = 3
         success = .false.
         do ii = 1,3
            jj = int(rand()*ntry) + 1
            ipairs = lst4(ip(itry(jj),:))
            RING_if: if (OH_groups_OK(ipairs(1),ipairs(2)).and.OH_groups_OK(ipairs(3),ipairs(4))) then
               ! Check if the 2 pairs of OH are are close enough
               r1 = rxyz(ipairs(1),:) - rxyz(ipairs(2),:)
               r2 = rxyz(ipairs(3),:) - rxyz(ipairs(4),:)
               call pbc(r1)
               call pbc(r2)
               rl1 = dot_product(r1,r1)
               rl2 = dot_product(r2,r2)
               OH_if: if ( max(rl1,rl2) < rmax2 ) then
                  ! Check if the 1st two Si are are close enough
                  iSi1 = 0
                  do i = 1,2
                     k = proximity(ipairs(1),i)
                     if (k == 0) cycle
                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                     iSi1 = k
                  end do
                  iSi2 = 0
                  do i = 1,2
                     k = proximity(ipairs(2),i)
                     if (k == 0) cycle
                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                     iSi2 = k
                  end do
                  r1 = rxyz(iSi1,:) - rxyz(iSi2,:)
                  call pbc(r1)
                  rl1 = dot_product(r1,r1)
                  Si_1_if: if ( rl1 < rmaxSiSi2 ) then
                     ! Check if the 2nd two Si are are close enough
                     iSi1 = 0
                     do i = 1,2
                        k = proximity(ipairs(3),i)
                        if (k == 0) cycle
                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                        iSi1 = k
                     end do
                     iSi2 = 0
                     do i = 1,2
                        k = proximity(ipairs(4),i)
                        if (k == 0) cycle
                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
                        iSi2 = k
                     end do
                     r2 = rxyz(iSi1,:) - rxyz(iSi2,:)
                     call pbc(r2)
                     rl2 = dot_product(r2,r2)
                     Si_2_if: if ( rl2 < rmaxSiSi2 ) then
                        success = .true.
                        exit
                     else Si_2_if
                        ! print *,'Si - Si 2 failed'
                     end if Si_2_if
                  else Si_1_if
                     ! print *,'Si - Si 1 failed'
                  end if Si_1_if
               else OH_if
                  ! print *,'OH constraint failed'
               end if OH_if
            else RING_if
               ! print *,'Ring constraint failed'
               itry(jj) = itry(ntry)
               ntry = ntry - 1
            end if RING_if
         end do
         if (.not.success) then
            call restore_config()
            call remove_from_list(lstSi,iat)
            GOTO 100
!            print *,'could not find 2 suitable pairs to condense'
!            cycle etching
         end if

         r1 = rxyz(ipairs(2),:) - rxyz(ipairs(1),:)
         call pbc(r1)
         rxyz(ipairs(1),:) = rxyz(ipairs(1),:) +  r1*0.5_wp
         call pbc(rxyz(ipairs(1),:))

         r2 = rxyz(ipairs(4),:) - rxyz(ipairs(3),:)
         call pbc(r2)
         rxyz(ipairs(3),:) = rxyz(ipairs(3),:) +  r2*0.5_wp
         call pbc(rxyz(ipairs(3),:))

         if (maxval(ipairs(1:2)) > maxval(ipairs(3:4))) then
            call condensation(ipairs(1),ipairs(2))
            call condensation(ipairs(3),ipairs(4))
         else
            call condensation(ipairs(3),ipairs(4))
            call condensation(ipairs(1),ipairs(2))
         end if


         call bond_list
         if (.not.allocated(bondl)) allocate(bondl(nbondtot))

         do ib = 1,nbondtot
            r1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
            call pbc(r1)
            bondl(ib) = sqrt(dot_product(r1,r1))
            ! print *,ib,ibond(1,ib),ibond(2,ib), bondl(ib)
            ! energy = 0.5_wp*KSiO*(r1sq - ASiO)**2
         end do
         ibmax = maxloc(bondl(1:nbondtot))
!         print *,bondl(ibmax(1))
         bondlmax = bondl(ibmax(1))

         call store_config1()
         call bond_switch_O(ipairs(1),success)
         call bond_switch_O(ipairs(3),success)

         call bond_angle_list()
         call bond_list
         if (.not.allocated(bondl)) allocate(bondl(nbondtot))

         do ib = 1,nbondtot
            r1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
            call pbc(r1)
            bondl(ib) = sqrt(dot_product(r1,r1))
         end do
         ibmax = maxloc(bondl(1:nbondtot))
         call pbc(r1)
         if (bondl(ibmax(1)) > bondlmax) then
            call restore_config1()
         else
            print *,'switch accepted'
         end if






!call check_proximity(ok,ii)
!if (.not.ok) then
!print *,'just after condensation 2'
!print *,'proximity array is inconsistent ',ii
!stop
!end if
!do i = 1,natom
!ia = atom(i)
!select case(ia)
!case(iSilicon)  ! Check Si is bonded to 4 atoms
!   do j = 1,ncmax(ia)
!      k = proximity(i,j)
!      if (k == 0) STOP 'ERROR 01'
!      if (atom(k) ==iSilicon) STOP 'ERROR 02'
!      if (atom(k) ==iHydrogen) STOP 'ERROR 03'
!   end do
!case(iOxygen)  ! Bridging Oxygen must be bonded to Si
!   i1 = proximity(i,1)
!   if (i1 == 0) then
!      print *,'i = ',i
!      print *,'ipairs = ',ipairs
!      STOP 'ERROR 04 B'
!   end if
!   a1 = atom(proximity(i,1))
!   if (a1 /= iSilicon) STOP 'ERROR 05'
!      i2 = proximity(i,2)  !  must be bonded to  2 Silicons
!      if (i2 == 0) STOP 'ERROR 06'
!      a2 = atom(proximity(i,2))
!      if (a2 /= iSilicon) STOP 'ERROR 07'
!   if (any(proximity(i,3:) /= 0)) STOP 'ERROR 08'
!case default
!   stop 'unknown atom type'
!end select
!end do


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

! check the clusters
         atomL = 0
         call Init_HKNonLattice(natom)
         call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
         !print *,'n_cluster = ',n_cluster
         if (n_cluster > 1) then
!           call remove_from_list(lstSi,iat)
!           if (atom(iat) ==iSilicon) then
!              where(lstSi%i == (natom+1)) lstSi%i = iat
!           end if
            call restore_config()
            cycle etching
         end if

         ! check for errors
         if (ANY(proximity(1:natom,:) > natom)) then
            print *,'error: proximity(1:natom,:) > natom '
            do i = 1,natom
               if (ANY(proximity(i,:) > natom)) then
               print *,i,' : ',proximity(i,:)
               end if
            end do
            stop
         end if
         if (ANY(proximity(natom + 1:,:) /= 0)) then
            print *,'error: proximity(natom + 1:,:) /= 0 '
            stop
         end if

         fdens_wt = frac_dens_wt()
         if (mod(isil,10) == 0) print '(i6,f9.5)',isil,fdens_wt

         if (fdens_wt <= frem) exit etching

      end do etching
!
!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      call bond_list
      proximity2 = proximity
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (dot_product(r3,r3) > boxl2*boxl2 ) then
            do j = 1,4
               if (proximity2(ibond(1,i),j) == ibond(2,i)) then
                   proximity2(ibond(1,i),j) = 0
               end if
            end do
            do j = 1,4
               if (proximity2(ibond(2,i),j) == ibond(1,i)) then
                   proximity2(ibond(2,i),j) = 0
               end if
            end do
         end if
         if (ABS(r3(1)) > boxl2) then
            n_crossbond_x = n_crossbond_x + 1
            crossbond_x(n_crossbond_x) = i
         end if
         if (ABS(r3(2)) > boxl2) then
            n_crossbond_y = n_crossbond_y + 1
            crossbond_y(n_crossbond_y) = i
         end if
         if (ABS(r3(3)) > boxl2) then
            n_crossbond_z = n_crossbond_z + 1
            crossbond_z(n_crossbond_z) = i
         end if
      end do
      atomL = 0

      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity2(1:natom,:),n_cluster,atomL)

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
!
! Finding a cluster which crosses all faces
!
       icx_n = 0
       icy_n = 0
       icz_n = 0
       do ibx = 1,n_crossbond_x
          icx = crossbond_x(ibx)
          if ( atomL(ibond(1,icx)) == atomL(ibond(2,icx)) ) then
             icx_n = icx_n + 1
             cluster_list_x(icx_n) = atomL(ibond(1,icx))
          end if
       end do

       do iby = 1,n_crossbond_y
          icy = crossbond_y(iby)
          if ( atomL(ibond(1,icy)) == atomL(ibond(2,icy)) ) then
             icy_n = icy_n + 1
             cluster_list_y(icy_n) = atomL(ibond(1,icy))
         end if
       end do

       do ibz = 1,n_crossbond_z
          icz = crossbond_z(ibz)
          if ( atomL(ibond(1,icz)) == atomL(ibond(2,icz)) ) then
             icz_n = icz_n + 1
             cluster_list_z(icz_n) = atomL(ibond(1,icz))
          end if
       end do

      spanning_cluster = .false.
      outer: do icl = 1, icx_n
      do icly = 1, icy_n
      do iclz = 1, icz_n
         if ((cluster_list_x(icl) == cluster_list_y(icly)) .AND. &
             (cluster_list_x(icl) == cluster_list_z(iclz))) then
             spanning_cluster = .true.
             exit outer
         end if
      end do
      end do
      end do outer

      do icly = 1,natom
         if ( atomL(icly) == cluster_list_x(icl) ) exit
      end do

      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      j = natom
      do
         if (atomL(j) /= atomL(icly)) then
print *,'deleting atoms not in 1st spanning cluster'
            call delete_atom(j)
            atomL(j) = atomL(icly)
         end if
         j = j - 1
         if (j == 0) exit
      end do

      do i = 1,natom
         if ((atom(i) == iSilicon).and.(atomL(i) == atomL(icly))) then
         do j = 1,4
            ii = proximity(i,j)
            if (ii == 0) STOP '5 error in proximity'
         end do
         end if
      end do

      write(*,*) 'natom ',natom
      write(*,*) 'irand_seed ',irand_seed
      if (spanning_cluster) then
         imve = 77
         write(unit=ctmp,fmt='(i6.6)') irand_seed0
         write(unit=c6,fmt='(i6.6)') natom
         write(unit=c5,fmt='(f5.3)') frem
         open(unit=imve,file=trim(outfile)//'_'//c5//'_'//trim(c6)//'_'//trim(ctmp)//'.pdb')
         call write_frame(imve,1,natom)
      else
         print *,'A spanning cluster was NOT generated'
      end if


! check everything is ok
! check each Si has 4 bonds
      do i = 1,natom
         if ((atom(i) == iSilicon)) then
         do j = 1,4
            ii = proximity(i,j)
            if (ii == 0) STOP '5 error in proximity'
         end do
         end if
      end do
      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
      call check_proximity2(ok,ii)
      if (.not.ok) then
         call write_proximity(6,1)
         print *,'proximity array is inconsistent (2) ',ii
         stop
      end if
!
! calculate the density
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      write(*,*) 'density % = ',100.0*frac_dens_wt()
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
      STOP 'normal termination'

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
         case default
            stop 'unknown atom type'
         end select
      end do
      frac_dens_wt = mass/(8000*mSi + 16000*mOx)
   END FUNCTION frac_dens_wt

   SUBROUTINE store_config()
      atom_old = atom
      rxyz_old = rxyz
      proximity_old = proximity
      natom_old = natom
   END SUBROUTINE

   SUBROUTINE restore_config()
      atom = atom_old
      rxyz = rxyz_old
      proximity = proximity_old
      natom = natom_old
   END SUBROUTINE

   SUBROUTINE store_config1()
      atom_old1 = atom
      rxyz_old1 = rxyz
      proximity_old1 = proximity
      natom_old1 = natom
   END SUBROUTINE

   SUBROUTINE restore_config1()
      atom = atom_old1
      rxyz = rxyz_old1
      proximity = proximity_old1
      natom = natom_old1
   END SUBROUTINE

END PROGRAM SIMBOX

