!!>include 'precision_mod.f90'

MODULE precision_mod
   implicit none
!  integer, parameter:: sp = kind(1.0)
!  integer, parameter:: dp = kind(1.0d0)
   integer, parameter:: sp = selected_real_kind(6,30)
   integer, parameter:: dp = selected_real_kind(15,300)
   integer, parameter:: qp_preferred = selected_real_kind(30,1000)
   integer, parameter:: qp = (1 + sign(1,qp_preferred))/2*qp_preferred+ &
                              (1 - sign(1,qp_preferred))/2*dp
!
   integer,parameter,public:: wp = dp
   integer,parameter,public:: i4b = selected_int_kind(9)
END MODULE precision_mod

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

!!>include 'global_vars.f90'

MODULE global_vars
   USE precision_mod, only: wp
   real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, etot, kBoltzT
   integer:: natom,nseed = 0,natom_max
   integer:: n_GAS,n_atom_GAS,n_Ne = 0,n_atom_Ne=0
   integer:: n_He = 0,n_atom_He=0
   integer:: nattached, nrelax
END MODULE

!!>include 'HKNonLattice2.f90'

MODULE HKNonLattice_mod
   implicit none
   integer:: n_cluster = 0,n_cluster_old = 0
   integer,parameter,private:: ncmx = 2000
   integer,allocatable:: atomL(:),atomL_old(:)
   integer,allocatable,private:: cluster(:,:)
   integer,private:: clusterC(ncmx)
CONTAINS

   SUBROUTINE Init_HKNonLattice(natom_max)
      integer,intent(in):: natom_max
      if (.not.allocated(atomL)) then
         allocate(atomL(natom_max),atomL_old(natom_max))
         !allocate(cluster(ncmx,natom_max))
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

! some routines used mainly for debugging
!   SUBROUTINE analyse_cluster(natom)
!      implicit none
!      integer,intent(in):: natom
!      integer:: i,L
!      clusterC(1:n_cluster) = 0
!      cluster(1:n_cluster,1:ncmx) = 0
!      do i = 1,natom
!         L = AtomL(i)
!         clusterC(L) = clusterC(L) + 1
!         cluster(L,clusterC(L)) = i
!      end do
!      ! some consistency checks
!      do i = 1,n_cluster
!         if (clusterC(i) /= count(atomL(1:natom) == i)) then
!            stop 'error in analyse_cluster'
!         end if
!      end do
!      if (sum(clusterC(1:n_cluster)) /= natom) then
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
!         if (clusterC(i) /= count(atomL(1:natom) == i)) then
!            stop 'error in print_cluster'
!         end if
!      end do
!   END SUBROUTINE print_cluster

END MODULE HKNonLattice_mod

!!>include 'seaton_mod.f90'

MODULE seaton_mod
   USE precision_mod
   USE global_vars, only: angstrom, K_ev, qstar, pi
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
   integer,parameter:: in_N2 = 8
   integer,parameter:: iO_O2 = 9
   integer,parameter:: iHe = 10
   integer,parameter:: iNe = 11
   integer,parameter:: ntyplj = 11  ! number of LJ types
   integer,parameter:: iN_N2charge = ntyplj+1
   integer,parameter:: iO_O2charge = ntyplj+2
   real(wp),parameter:: eps_O_sil   = 185.0_wp*K_ev, sig_O_sil  = 2.708_wp/angstrom,q_O_sil = -0.64025_wp*qstar
   real(wp),parameter:: eps_Si_sil  =  0.0_wp*K_ev, sig_Si_sil = 0.0_wp/angstrom,  q_Si_sil = -2.0_wp*q_O_sil
   real(wp),parameter:: eps_O_OH_sil = 185.0_wp*K_ev, sig_O_OH_sil = 3.0_wp/angstrom, q_O_OH_sil = -0.533_wp*qstar
   real(wp),parameter:: eps_H_OH_sil =  0.0_wp*K_ev, sig_H_OH_sil = 0.0_wp/angstrom, q_H_OH_sil =  0.206_wp*qstar
   real(wp),parameter:: eps_H_H2O   =  0.0_wp*K_ev,  sig_H_H2O = 0.0_wp/angstrom,   q_H_H2O =  0.417_wp*qstar
   real(wp),parameter:: eps_O_H2O   = 76.58_wp*K_ev, sig_O_H2O = 3.1506_wp/angstrom,q_O_H2O = -2.0_wp*q_H_H2O
   real(wp),parameter:: eps_O_CO2   = 82.997_wp*K_ev,sig_O_CO2 = 3.064_wp/angstrom, q_O_CO2 = -0.33225_wp*qstar
   real(wp),parameter:: eps_C_CO2   = 29.999_wp*K_ev,sig_C_CO2 = 2.785_wp/angstrom, q_C_CO2 = -2.0_wp*q_O_CO2
   real(wp),parameter:: eps_n_N2    = 34.897_wp*K_ev, sig_n_N2 = 3.3211_wp/angstrom, q_n_N2 = -0.5475_wp*qstar
   real(wp),parameter:: eps_O_O2    = 43.183_wp*K_ev, sig_O_O2 = 3.1062_wp/angstrom, q_O_O2 = -0.3577_wp*qstar
   real(wp),parameter:: eps_He      = 11.0_wp*K_ev,   sig_He   = 2.64_wp/angstrom,     q_He =  0.0_wp*qstar
   real(wp),parameter:: eps_Ne      = 50.0_wp*K_ev,   sig_Ne   = 2.76_wp/angstrom,     q_Ne =  0.0_wp*qstar
   real(wp),parameter:: eps_He_Si   = 39.0_wp*K_ev,  sig_He_Si = 2.39_wp/angstrom
   real(wp),parameter:: eps_He_O    = 29.0_wp*K_ev,  sig_He_O  = 2.83_wp/angstrom
   real(wp),parameter:: eps_Ne_Si   = 84.0_wp*K_ev,  sig_Ne_Si = 3.45_wp/angstrom
   real(wp),parameter:: eps_Ne_O    = 54.0_wp*K_ev,  sig_Ne_O  = 2.92_wp/angstrom
!  bond lengths
   real(wp),parameter:: bondl_SiO = 1.60_wp/angstrom
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
   real(wp),parameter:: theta_CO2 = pi
   real(wp),parameter:: K_CO2 = 5.331    ! eV
   real(wp),parameter:: KOO = 366.142_wp    ! eV nm^-2
   real(wp),parameter:: KNN = 366.142_wp    ! eV nm^-2
   real(wp),parameter:: KCO2 = 366.142_wp   ! eV nm^-2
   real(wp),parameter:: AOO = bondl_O2     ! nm
   real(wp),parameter:: ANN = bondl_N2     ! nm
   real(wp),parameter:: ACO2 = bondl_CO2    ! nm
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
                        eps_n_N2, &
                        eps_O_O2, &
                        eps_He, &
                        eps_Ne /)
   real(wp),parameter:: sigi(0:ntyplj) = (/ &
                        sig_O_sil, &
                        sig_Si_sil, &
                        sig_O_OH_sil, &
                        sig_H_OH_sil, &
                        sig_O_H2O, &
                        sig_H_H2O, &
                        sig_O_CO2, &
                        sig_C_CO2, &
                        sig_n_N2, &
                        sig_O_O2, &
                        sig_He, &
                        sig_Ne /)
   real(wp),parameter:: sigi2(0:ntyplj) = sigi*0.5_wp
! Charges
   real(wp),parameter:: q_seaton(0:ntyplj+2) = (/ &
                        q_O_sil, &
                        q_Si_sil, &
                        q_O_OH_sil, &
                        q_H_OH_sil, &
                        q_O_H2O, &
                        q_H_H2O, &
                        q_O_CO2, &
                        q_C_CO2, &
                        q_n_N2, &
                        q_O_O2, &
                        q_He, &
                        q_Ne, &
                -2.0_wp*q_n_N2, &
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
   real(wp):: energy
!  united atom H2O
   real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
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
   real(wp),parameter:: sigma(0:5) = (/ sig_O_sil, sig_Si_sil,sigma_OH, sig_H_OH_sil,sig_O_H2O, sig_H_H2O /)
   real(wp),parameter:: sigma_2(0:5) = sigma*0.5_wp
END MODULE

!!>include 'atom_types.f90'

MODULE atom_types
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
   real(wp),parameter:: mSi = 28.06_wp
   real(wp),parameter:: mOx = 16.00_wp
   real(wp),parameter:: mOH = 17.00_wp
   real(wp),parameter:: mH =  1.00_wp
   real(wp),parameter:: mN = 14.00_wp
   real(wp),parameter:: mC = 12.0107_wp
   real(wp),parameter:: mHe = 4.002602_wp
   real(wp),parameter:: mNe = 20.1797_wp
   real(wp),parameter:: amass(0:ntyp) = (/ mOx,mSi,mOH,mH,mOx,mH /)
   integer,allocatable:: atom(:),He_atom(:),Ne_atom(:)
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

!!>include 'coordinates_mod.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:)
   real(wp),allocatable:: He_xyz(:,:),r_He_uc(:,:),dr_He(:,:)
   real(wp),allocatable:: He_fxyz(:,:),He_vxyz(:,:)
   real(wp),allocatable:: Ne_xyz(:,:),r_Ne_uc(:,:),dr_Ne(:,:)
   real(wp),allocatable:: Ne_fxyz(:,:),Ne_vxyz(:,:)
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
      USE atom_types, only: atom_name,atom
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
            if ( any(proximity(i,3:) /= 0) ) GOTO 100
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

!!>include 'Keating_parameters_vonAlfthan.f90'

MODULE Keating_parameters
!
! Parameters for the simplified Keating potential.
! Source: S. von alfthan et al., Phys. Rev. B 68, 073203, (2003)
!
   USE precision_mod, only: wp
   implicit none
   real(wp),parameter:: KSiSi = 908.0_wp    ! eV nm^-2
   real(wp),parameter:: KSiO = 2700.0_wp
   real(wp),parameter:: ASiSi = 0.235_wp    ! nm
   real(wp),parameter:: ASiO = 0.161_wp
   real(wp),parameter:: KSiSiSi = 3.58_wp   ! eV
   real(wp),parameter:: KSiSiO = 3.93263_wp
   real(wp),parameter:: KOSiO = 4.32_wp
   real(wp),parameter:: KSiOSi = 2.0_wp
   real(wp),parameter:: kbond(1:3) = (/ KSiO,KSiSi,KSiO /)
   real(wp),parameter:: abond(1:3) = (/ ASiO,ASiSi,ASiO /)
   real(wp),parameter:: acosang(0:2) = (/ 1.0_wp,1.0_wp/3.0_wp,1.0_wp /)
END MODULE

!!>include 'bond_list_mod.f90'

MODULE bond_list_mod
   USE precision_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   integer:: nbondtot,nang
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types
      USE Keating_parameters
      integer:: ia1,ia2,j
!
      nbondtot = 0
      if (.not.allocated(ibond)) then
         allocate( ibond(2,natom_max*4) )
      end if
      if (.not.allocated(nbond)) then
         allocate( nbond(natom_max) )
      end if
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
            end if
         end do
      end do
   END SUBROUTINE bond_list

   SUBROUTINE set_bond_angle_lists()
      USE global_vars, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types
      USE Keating_parameters
      integer:: i,ia1,ia2,ia3,ib,ic
!
      if (.not.allocated(ibond)) then
         allocate( ibond(2,natom_max*4) )
         allocate( iang(3,natom_max*6) )
      end if
      allocate(kang(6*natom_max),ctheta(6*natom_max))
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
            end if
         end do
      end do
      end do
      do i = 1, nang
         if (atom(iang(2,i)) == iSilicon) then
            kang(i) = KOSiO
            ctheta(i) = -1.0_wp/3.0_wp
         else
            kang(i) = KSiOSi
            ctheta(i) = -1.0_wp
         end if
      end do
   END SUBROUTINE set_bond_angle_lists

END MODULE bond_list_mod

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
      USE global_vars
      USE atom_types
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + n_He + n_Ne + 2
            write(iu,'(a,6f14.8)')'Gases+Silica ',boxl*Angstrom
         else
            write(iu,*) n_He + n_Ne + 2
            write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
         end if
      else
         write(iu,*) n_He + n_Ne + 2
         write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
      end if
      do i = 1,n_He
         write(iu,110)'He ',He_xyz(i,:)*Angstrom
      end do
      do i = 1,n_Ne
         write(iu,110)'Ne ',Ne_xyz(i,:)*Angstrom
      end do
!
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,110) atom_name(atom(i)),(rxyz(i,j)*Angstrom,j = 1,NDIM)
         end do
         end if
      end if
      write(iu,110)'Ar ',  boxl2*Angstrom
      write(iu,110)'Ar ', -boxl2*Angstrom
      write(iu,*)
      close(iu,status='keep')
110   format(a,6(1x,f12.6))
   END SUBROUTINE write_box_frame

   SUBROUTINE write_vel(iu,nmol,print_all)
      USE global_vars
      USE atom_types
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j,k
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + n_He +n_Ne
            write(iu,'(a,6f14.8)')'Gases+Silica ',boxl*Angstrom
         else
            write(iu,*) n_He +n_Ne
            write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
         end if
      else
         write(iu,*) nmol + nmol
         write(iu,'(a,6f14.8)')'Gases ',boxl*Angstrom
      end if
      do i = 1,n_He
         write(iu,*) (He_vxyz(i,k),k = 1,NDIM)
      end do
      do i = 1,n_Ne
         write(iu,*) (Ne_vxyz(i,k),k = 1,NDIM)
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
      USE global_vars
      USE atom_types
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      integer:: i,j,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp
      do i = 1,n_He
         read(iu,*) He_vxyz(i,:)
      end do
      do i = 1,n_Ne
         read(iu,*) Ne_vxyz(i,:)
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
      USE global_vars
      USE atom_types
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      real(wp):: r(NDIM)
      integer:: i,j,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp
      do i = 1,n_He
         read(iu,*) ctmp,r
         He_xyz(i,:) = r/Angstrom
      end do
      do i = 1,n_Ne
         read(iu,*) ctmp,r
         Ne_xyz(i,:) = r/Angstrom
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

END MODULE BOX_FRAMES_MOD

!!>include 'comvel_mod.f90'

MODULE comvel_mod
   USE precision_mod
   USE global_vars
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

!!>include 'force_keating.f90'

MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none
CONTAINS

   SUBROUTINE force_keating
      USE bond_list_mod
      USE Keating_parameters
      USE coordinates_mod
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: k,dist,rl1(NDIM),r1sq,ktheta,rl2(NDIM),r2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
      fxyz = 0.0_wp
      energy = 0.0_wp

!     Force routine for Keating bond stretching
      do ib = 1,nbondtot
         k = KSiO
         dist = ASiO
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) - k*(r1sq - dist)*rl1(1:NDIM)/r1sq
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) + k*(r1sq - dist)*rl1(1:NDIM)/r1sq
         energy = energy + 0.5_wp*k*(r1sq - dist)**2
      end do

!     Force routine for Keating bond angle bending
      do ia = 1,nang
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:NDIM) - rxyz(ia2,1:NDIM)
         rl2 = rxyz(ia3,1:NDIM) - rxyz(ia2,1:NDIM)
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
         fxyz(ia2,1:NDIM) = fxyz(ia2,1:NDIM) - coef*(-(rl1 + rl2)/coef2 &
                                       + rl1l2*rl1/coef1 &
                                       + rl1l2*rl2/coef3)
         fxyz(ia1,1:NDIM) = fxyz(ia1,1:NDIM) - coef*(rl2/coef2 - rl1l2*rl1/coef1)
         fxyz(ia3,1:NDIM) = fxyz(ia3,1:NDIM) - coef*(rl1/coef2 - rl1l2*rl2/coef3)
         energy = energy + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE force_keating

END MODULE force_keating_mod

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
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: r(:,:),v(:,:)
      real(wp),intent(inout):: dr(:,:)
      real(wp):: dt2,dtsq2
      dt2 = dt*0.5_wp
      dtsq2 = dt*dt2
      dr = dt*v + dtsq2*f*massi
      r = r + dr
      v = v + dt2*f*massi
   END SUBROUTINE gas_vverlet_a

   SUBROUTINE gas_vverlet_b(dt,massi,v,f)
      real(wp),intent(in):: dt,massi(:,:),f(:,:)
      real(wp),intent(inout):: v(:,:)
      v = v + dt*0.5_wp*f*massi
   END SUBROUTINE gas_vverlet_b

END MODULE vverlet_mod

!!>include 'charges2_mod.f90'

MODULE charges_mod
   USE precision_mod
   USE global_vars, only: nseed
   USE atom_types
   USE seaton_mod, only: qi => q_seaton
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
      USE atom_types
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


!!>include 'nlist_mod.f90'

MODULE NLIST_MOD
   USE precision_mod, only: wp
   USE global_vars, only: natom_max,natom
   USE coordinates_mod
   implicit none
   PRIVATE
   integer,public:: neigh
   real(wp),private:: delx,dely,delz,delmin
   real(wp),private:: delxi,delyi,delzi
   integer,public:: ncelx,ncely,ncelz,ncelt
   integer,parameter:: neighmx = 350
   integer,allocatable,public:: hoc(:),hoc_old(:)
   integer,public:: ncell(neighmx)
   integer,allocatable,public:: ll_old(:),ll(:),lr(:)
   public:: NEIGCELL,cell,INIT_NLIST,NEW_NLIST,print_cell
   public:: push_cell,pop_cell,Z_NEIGH_CELL,nlayers_ll
!
!  ll(i)      : linked list particle i
!  hoc(ic)    : head of chain cell ic
!  ncelx,y,z  : number of cells in x, y or z direction
!  ncelt      : total number of cells
!  neigh      : number of cells for interactions
!
CONTAINS

   PURE FUNCTION nlayers_ll(r)
      real(wp),intent(in):: r
      integer:: nlayers_ll
      nlayers_ll = int(r/delmin) + 1
   END FUNCTION

   SUBROUTINE INIT_NLIST(Lx,Ly,Lz,Rc)
      real(wp),intent(in):: Lx,Ly,Lz,Rc
      allocate(ll_old(natom_max),ll(natom_max),lr(natom_max))
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

   SUBROUTINE NEW_NLIST
!     makes a new neighbour list using the linked-list algorithm
      integer:: i,ic
      HOC(1:NCELT) = 0  ! initialize the head-of-chain
      ! make linked list:
      do i = 1,natom
         ! determine cell number
         ic = CELL(rxyz(i,:))
         if (ic < 1 .or. ic > ncelt) stop 'NEW_NLIST: cell number out of bounds'
         ! update linked - list and head of chain
         LL(i) = HOC(ic)
         if (HOC(ic) /= 0) LR(HOC(ic)) = i
         HOC(ic) = i
      end do
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


!!>include 'repul_force_NL_vonAlfthan.f90'

MODULE repul_energy_mod
   USE precision_mod, only: wp
   implicit none
   real(wp):: repulenergy
CONTAINS

   SUBROUTINE repul_energy(ifirst,ilast)
      USE coordinates_mod
      USE atom_types, only: atom,iHydrogen
      USE connectivity_mod, only: nearest_neighbor2
      USE NLIST_MOD
      integer,intent(in):: ifirst,ilast
      real(wp):: rr(NDIM),r2,r(NDIM)
      integer:: M,L,ic,nc,inn,nn,neigbors(27)
      real(wp),parameter:: rden = 0.26_wp  ! nm
      real(wp),parameter:: rden2 = rden**2
      real(wp),parameter:: ff = 8000.0_wp  ! eV nm^-4
!
      repulenergy = 0.0_wp
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
               repulenergy = repulenergy + (r2 - rden2)**2
               fxyz(L,:) = fxyz(L,:) + 2.0_wp*ff*(r2 - rden2)*(rr)
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
      repulenergy = repulenergy*0.5_wp*ff
   END SUBROUTINE

   SUBROUTINE repul_energy2(ifirst,ilast)
      USE coordinates_mod
      USE atom_types, only: atom,iHydrogen
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
               repulenergy = repulenergy + 0.5_wp*ff*(r2 - rden2)**2
               fxyz(L,1) = fxyz(L,1) + 2.0_wp*ff*(r2 - rden2)*dx
               fxyz(L,2) = fxyz(L,2) + 2.0_wp*ff*(r2 - rden2)*dy
               fxyz(L,3) = fxyz(L,3) + 2.0_wp*ff*(r2 - rden2)*dz
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
   END SUBROUTINE
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

!!>include 'insert_mod.f90'

MODULE insert_mod
    USE precision_mod, only: wp
    USE global_vars
    USE nlist_mod
    USE atom_types
    USE coordinates_mod
    USE seaton_mod
    USE rand_mod
    USE quat2mat_mod
    USE ran_point_sphere
    implicit none

CONTAINS

   SUBROUTINE insert_mol_GAS(n_mol,real_xyz,ntrial,rad1,fover0)
      integer,intent(in):: ntrial
      integer,intent(inout):: n_mol
      real(wp),intent(in):: rad1
      real(wp),intent(out):: real_xyz(:,:)
      real(wp),intent(in),optional:: fover0
      real(wp):: x,y,z,fover
      real(wp):: mol_vec1(NDIM)
      integer:: it,im,n_mol_requested
!
      if (present(fover0)) then
         fover = fover0
      else
         fover = 1.0_wp
      end if
!
      n_mol_requested = n_mol
      n_mol = 0
      mol_loop: do im = 1, n_mol_requested
      insert_loop: do it = 1,ntrial

         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2(1)
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2(2)
         z = ( 1.0_wp - 2.0_wp*rand() )*boxl2(3)
         mol_vec1 = (/ x,y,z /)

! first check if there is a large overlap
         if (.not.overlap(fover,rad1,mol_vec1)) then
            n_mol = n_mol + 1
            real_xyz(n_mol,:) = mol_vec1
            write(*,*) im,n_mol
            cycle mol_loop
         end if

      end do insert_loop
      end do mol_loop
   END SUBROUTINE insert_mol_GAS


   FUNCTION overlap(fr,rad,mol_vec1)
      logical:: overlap
      real(wp),intent(in):: fr,rad,mol_vec1(:)
      integer:: jj,j,ic,nc
      real(wp):: dr(NDIM),sigj2
      overlap = .false.
      ic = CELL( mol_vec1(1:NDIM) )
      call NEIGCELL(ic,2,neigh,ncell)
      cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         j = HOC(nc)
         cell_atom_loop: do while (j /= 0)
            sigj2 = sigi2(atom(j))
            if (atom(j) == iSilicon) then
               sigj2 = 2.07_wp/angstrom
            end if
            dr(1:NDIM) = mol_vec1(1:NDIM) - rxyz(j,1:NDIM)
            call pbc(dr)
            if (dot_product(dr,dr) < (fr*(rad + sigj2))**2) then
               overlap = .true.
               RETURN
            end if
            j = LL(j)
         end do cell_atom_loop
      end do cell_loop
   END FUNCTION


   FUNCTION overlap2(fr,rad,mol_vec1)
      logical:: overlap2
      real(wp),intent(in):: fr,rad,mol_vec1(:)
      integer:: j
      real(wp):: dr(NDIM),sigj2
      overlap2 = .false.
      do j = 1,natom
         sigj2 = sigi2(atom(j))
         if (atom(j) == iSilicon) then
            sigj2 = 2.07_wp/angstrom
         end if
         dr(1:NDIM) = mol_vec1(1:NDIM) - rxyz(j,1:NDIM)
         call pbc(dr)
         if (dot_product(dr,dr) < (fr*(rad + sigj2))**2) then
            overlap2 = .true.
            RETURN
         end if
      end do
   END FUNCTION

END MODULE insert_mod

!!>include 'lj_el_mod_O2_N2_CO2.f90'

MODULE Lj_el_mod
    USE precision_mod
    USE global_vars
    USE coordinates_mod
    USE Keating_parameters
    USE charges_mod
    USE atom_types
    USE seaton_mod
    USE constants_mod
    USE nlist_mod
    implicit none
    save
    real(wp),allocatable:: aij(:,:),bij(:,:)
    real(wp),allocatable:: vdwcut(:,:),vdwsf(:,:),vdwcut_mod(:,:),vdwsf_mod(:,:)
    real(wp),parameter:: rcut = 12.0_wp/angstrom
    real(wp),parameter:: rcut2 = rcut**2
    real(wp),parameter:: rcut6 = 1.0_wp/(rcut2*rcut2*rcut2)
    real(wp),parameter:: rcut7 = rcut6/rcut
    real(wp),parameter:: rcut13 = rcut7*rcut7*rcut

CONTAINS

   SUBROUTINE LJ_INIT
      real(wp):: eij,rstij
      integer:: i,j
      allocate( aij(0:ntyplj,0:ntyplj) )
      allocate( bij(0:ntyplj,0:ntyplj) )
      write(*,*) 'ntyplj = ',ntyplj
      allocate( vdwcut(0:ntyplj,0:ntyplj) )
      allocate( vdwsf(0:ntyplj,0:ntyplj) )
      allocate( vdwcut_mod(0:ntyplj,0:ntyplj) )
      allocate( vdwsf_mod(0:ntyplj,0:ntyplj) )
!
print *,'Lennard-Jones'
      do j = 0,ntyplj    ! lj aij,bij parameters & shifted terms
      do i = 0,ntyplj
         eij = sqrt(epsi(i)*epsi(j))
         rstij = (sigi(i) + sigi(j))*0.5
if((j==iSilicon .and. i==iHe).or.(i==iSilicon .and. j==iHe))then
   eij = eps_He_Si
   rstij = sig_He_Si
print *,i,j,eij,rstij
end if
if((j==iSilicon .and. i==iNe).or.(i==iSilicon .and. j==iNe))then
   eij = eps_Ne_Si
   rstij = sig_Ne_Si
print *,i,j,eij,rstij
end if
if((j==iOxygen .and. i==iHe).or.(i==iOxygen .and. j==iHe))then
   eij = eps_He_O
   rstij = sig_He_O
print *,i,j,eij,rstij
end if
if((j==iOxygen .and. i==iNe).or.(i==iOxygen .and. j==iNe))then
   eij = eps_Ne_O
   rstij = sig_Ne_O
print *,i,j,eij,rstij
end if
if((j==iOxygenH .and. i==iHe).or.(i==iOxygenH .and. j==iHe))then
   eij = eps_He_O
   rstij = sig_He_O
end if
if((j==iOxygenH .and. i==iNe).or.(i==iOxygenh .and. j==iNe))then
   eij = eps_Ne_O
   rstij = sig_Ne_O
end if
         aij(i,j) = 4.0*(rstij**12)*eij
         bij(i,j) = 4.0*(rstij**6)*eij
         vdwcut(i,j) = rcut6*(aij(i,j)*rcut6 - bij(i,j))
         vdwsf(i,j) = rcut6*(-12.0_wp*aij(i,j)*rcut6 + 6.0*bij(i,j))/rcut
!        vdwcut_mod(i,j) = rcut_1_6*(aij(i,j)*rcut_1_6-bij(i,j))
!        vdwsf_mod(i,j) =  rcut_1_6*(-12.0_wp*aij(i,j)*rcut_1_6+6.0*bij(i,j))/rcut_1
!        write(*,*) i,j,aij(i,j),bij(i,j),vdwcut(i,j),vdwsf(i,j)
      end do
      end do
   END SUBROUTINE LJ_INIT

   SUBROUTINE FORCE_ENERGY_LJ_NLIST(gasfirst,gaslast,atom_gas,rxyz_gas,fxyz_gas,ULJEL)

      integer,intent(in):: gasfirst,gaslast
      integer,intent(in):: atom_gas(:)
      real(wp),intent(in):: rxyz_gas(:,:)
      real(wp),intent(inout):: fxyz_gas(:,:)
      real(wp),intent(out):: ULJEL(:)
      real(wp):: rr(NDIM),df(NDIM),r2,r1,rr6,dele,rr7,rr13
      integer,parameter:: nlm = 5
      integer:: L,M,ic,jj,nc,ncell((2*nlm + 1)**3),neigh,nl
      nl = nlayers_ll(rcut - tiny(1.0))
      if (nl > nlm) STOP 'nl > nlm '
      ULJEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
         ic = CELL( rxyz_gas(M,:) )
if (ic < 1) stop 'FORCE_ENERGY_LJ_NLIST: ic < 1'
if (ic > ncelt) stop 'FORCE_ENERGY_LJ_NLIST: ic > ncelt'
         call NEIGCELL(ic,nl,neigh,ncell)
         cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         L = HOC(nc)
         cell_atom_loop: do while (L /= 0)
            rr(1:NDIM) = rxyz_gas(M,:) - rxyz(L,:)
            call pbc(rr)
            r2 = dot_product(rr,rr)
            if (r2 < rcut2) then
            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7
            dele = rr6*(aij(atom(L),atom_gas(M))*rr6 - bij(atom(L),atom_gas(M))) &
                 -   vdwcut(atom(L),atom_gas(M)) - vdwsf(atom(L),atom_gas(M))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),atom_gas(M))*(rr7 - rcut7) &
               + 12.0_wp*aij(atom(L),atom_gas(M))*(rcut13 - rr13))*rr(1:NDIM)/r1
            fxyz(L,1:NDIM) = fxyz(L,1:NDIM) + df
            fxyz_gas(M,1:NDIM) = fxyz_gas(M,1:NDIM) - df
            ULJEL(M) = ULJEL(M) + dele
            end if
         L = LL(L)
         end do cell_atom_loop
         end do cell_loop
      end do mol_loop
      energy = energy + sum(ULJEL(gasfirst:gaslast))
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
      energy = energy + UEL
   END SUBROUTINE

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
!
      do i = 1,T%nmol
         T%hist(i,nt,:) = xyz(i,:)
      end do
!
      if (T%icount <= T%ntau) RETURN
      T%isamp = T%isamp + 1
      do L = 1,T%ntau
         m = nt - L + 1
         if (m <= 0) m = m + T%ntau
         do i = 1,T%nmol
            T%tcfmol(i,L) = T%tcfmol(i,L) &
                          + ( T%hist(i,nt,1) - T%hist(i,m,1) )**2 &
                          + ( T%hist(i,nt,2) - T%hist(i,m,2) )**2 &
                          + ( T%hist(i,nt,3) - T%hist(i,m,3) )**2
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
   END SUBROUTINE tetra_coords5

END MODULE utility_mod


PROGRAM SIMBOX
      USE precision_mod
      USE commandline_mod
      USE command_arg_mod
      USE precision_mod
      USE utility_mod
      USE files_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE keating_parameters
      USE nlist_mod
      USE repul_energy_mod
      USE charges_mod
      USE insert_mod
      USE tcf_mod
      USE msd_mod
      USE vcf_mod
      implicit none
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      integer:: number_moves, nequib, nsamp, nprint, nprint_all, nfreq_movie, n_short
      logical:: make_movie
      integer:: ntau_msd, ntau_vcf
      real(wp):: del_nlist, fover
      real(wp):: dt,qj,sumqj,T_kelvin,T_kinetic
      real(wp):: EK,Ek_He,Ek_Ne,Ek_tot
      real(wp):: dtreal,time,rr(NDIM),mtmp,ULJAtt
      real(wp),allocatable:: ULJEL(:)
      real(wp),allocatable:: massi(:,:),mass(:,:)
      real(wp),allocatable:: He_massi(:,:),He_mass(:,:)
      real(wp),allocatable:: Ne_massi(:,:),Ne_mass(:,:)
      real(wp):: rad_Ox(1),rad_Ne(1)
      integer:: n_OH
      integer:: j,i,istep,stat,nc
      integer:: itmp,narg,Ndf
      real(wp),allocatable:: coprxyz(:,:)
      real:: tp(0:10),sumtp(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5,c32
      character(len=132):: coordfile,restart_coordfile,restart_velfile,carg,ftyp*3,ctrlfile
      integer:: io_temp, io_energy, io_config, io_config_final, io_topol, io_ctrl
      integer:: io_restart_xyz, io_restart_v
      logical:: restart
      type(tcf_type):: msd_He, vcf_He, msd_Ne, vcf_Ne
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
      read(io_ctrl,*) c32, nprint              ; write(*,*) trim(c32), nprint
      read(io_ctrl,*) c32, nprint_all          ; write(*,*) trim(c32), nprint_all
      read(io_ctrl,*) c32, make_movie          ; write(*,*) trim(c32), make_movie
      read(io_ctrl,*) c32, nfreq_movie         ; write(*,*) trim(c32), nfreq_movie
      read(io_ctrl,*) c32, nsamp               ; write(*,*) trim(c32), nsamp
      read(io_ctrl,*) c32, ntau_msd            ; write(*,*) trim(c32), ntau_msd
      read(io_ctrl,*) c32, ntau_vcf            ; write(*,*) trim(c32), ntau_vcf
      read(io_ctrl,*) c32, N_He      ; write(*,*) trim(c32), n_He
      read(io_ctrl,*) c32, n_Ne      ; write(*,*) trim(c32), n_Ne
      read(io_ctrl,*) c32, del_nlist           ; write(*,*) trim(c32), del_nlist
      read(io_ctrl,*) c32, fover               ; write(*,*) trim(c32), fover
      del_nlist = del_nlist/Angstrom

      call next_command_argument(coordfile, "input file",stat)
      nc = len_trim(coordfile)
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
!
      allocate(He_xyz(1:n_He,1:NDIM))
      allocate(r_He_uc(n_He,3))
      allocate(dr_He(n_He,3))
      allocate(He_fxyz(1:n_He,1:NDIM))
      allocate(He_vxyz(1:n_He,1:NDIM))
      allocate(He_atom(1:n_He))
      allocate(He_mass(1:n_He,3))
      allocate(He_massi(1:n_He,3))
!
      allocate(Ne_xyz(1:n_Ne,1:NDIM))
      allocate(r_Ne_uc(n_Ne,3))
      allocate(dr_Ne(n_Ne,3))
      allocate(Ne_fxyz(1:n_Ne,1:NDIM))
      allocate(Ne_vxyz(1:n_Ne,1:NDIM))
      allocate(Ne_atom(1:n_Ne))
      allocate(Ne_mass(1:n_Ne,3))
      allocate(Ne_massi(1:n_Ne,3))
!

      allocate( ULJEL(max(n_He,n_Ne)) )

      select case(ftyp)
      case('xyz')
         read (io_config,*) natom
         read (io_config,*) boxl
         io_topol = 15
         call get_free_file_unit(io_topol)
         open(unit=io_topol,file=trim(coordfile(1:nc - 3))//'con')
      case('pdb')
         read (io_config,*) c6,natom  ;print *,'natom = ',natom
         read (io_config,*) c6,boxl   ;print '(a,6f14.8)',' boxl = ',boxl
         ! io_topol = io_config
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select

!
      natom_max = natom
      boxl = boxl/Angstrom
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!rcutel = boxl2
!rcutel2 = boxl2**2
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(fxyz(1:natom_max,3),vxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
      proximity = 0
      allocate(mass(natom_max,3),massi(natom_max,3))


      select case(ftyp)
      case('xyz')
         do i = 1,natom
            read (io_config,*) ctmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
         close(io_config)
      case('pdb')
         do i = 1,natom
            read (io_config,'(a6)',advance = 'no') c6
            read (io_config,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
      end select
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom

      do i = 1,natom
         if (atom(i) == iSilicon) mass(i,:) = mSi
         if (atom(i) == iOxygen)  mass(i,:) = mOx
         if (atom(i) == iOxygenH)  mass(i,:) = mOH
      end do
      massi = 1.0_wp/mass


      select case(ftyp)
      case('xyz')
         read(io_topol,*) itmp
         do i = 1,natom
            read(io_topol,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*j + 1:6 + 5*j + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(io_topol)
      case('pdb')
         do i = 1,natom
            read(io_config,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*j + 1:6 + 5*j + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(io_config)
      end select

      call assign_charge(1,natom)

      write(*,*) 'sum charge = ',sum(charge(1:natom))
      write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
      write(*,*) 'n Si = ',count(atom == iSilicon)
      write(*,*) 'n O  = ',count(atom == iOxygen)
      write(*,*) 'n OH = ',count(atom == iOxygenH)
      write(*,*) 'ntot = ',count(atom == iSilicon) + count(atom == iOxygen) + count(atom == iOxygenH)
!
      call INIT_NLIST(boxl(1),boxl(2),boxl(3),del_nlist)
      call NEW_NLIST

!
      do i = 1,n_He
        He_atom(i) = iHe
        He_mass(i,:) = mHe
      end do
      He_massi = 1.0_wp/He_mass
      do i = 1,n_Ne
        Ne_atom(i) = iNe
        Ne_mass(i,:) = mNe
      end do
      Ne_massi = 1.0_wp/Ne_mass

! Insertion

      if (N_He > 0) then
         call insert_mol_GAS(N_He,He_xyz,1000000,0.5_wp*Sig_He,fover); print *,'n_He = ',n_He
      end if
      if (n_Ne > 0) then
         call insert_mol_GAS(n_Ne,Ne_xyz,1000000,0.5_wp*Sig_Ne,fover); print *,'n_Ne = ',n_Ne
      end if
!
      if (restart) then
         call read_box_frame(io_restart_xyz,nmol = n_Ne,print_all = .true.)
         call read_vel(io_restart_v,nmol = n_Ne,print_all = .true.)
      end if
!
! Set up the uncorrected coords
      r_He_uc(1:n_He,:) = He_xyz(1:n_He,:)
      r_Ne_uc(1:n_Ne,:) = Ne_xyz(1:n_Ne,:)


      if (make_movie) then ! write out first frame of movie
         iframe = 0
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_He,print_all = .true.)
      end if
!========================================================================
!
      call LJ_INIT

      call comvel(natom,T_kelvin,vxyz,mass)
      call comvel(n_He,T_kelvin,He_vxyz,He_mass)
      call comvel(n_Ne,T_kelvin,Ne_vxyz,Ne_mass)
!
      print '(a/,3g18.8)','net Ne momentum', &
         sum(Ne_vxyz(:,1)*Ne_mass(:,1)),&
         sum(Ne_vxyz(:,2)*Ne_mass(:,2)),&
         sum(Ne_vxyz(:,3)*Ne_mass(:,3))

      print '(a/,3g18.8)','net solid momentum', &
         sum(vxyz(:,1)*mass(:,1)),&
         sum(vxyz(:,2)*mass(:,2)),&
         sum(vxyz(:,3)*mass(:,3))
      print '(a/,3g18.8)','net solid momentum', &
         dot_product(vxyz(:,1),mass(:,1)),&
         dot_product(vxyz(:,2),mass(:,2)),&
         dot_product(vxyz(:,3),mass(:,3))

! Rescale each gas to the set temperature
      Ek = kinetic_energy(natom,vxyz,mass)
      Ek_He = kinetic_energy(n_He,He_vxyz,He_mass)
      Ek_Ne = kinetic_energy(n_Ne,Ne_vxyz,Ne_mass)

      if (n_He > 0) then
         Ndf = 3*n_He
         T_kinetic = Ek_He*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [He] = ',T_kinetic
         He_vxyz = He_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_He = kinetic_energy(n_He,He_vxyz,He_mass)
         T_kinetic = Ek_He*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [He] = ',T_kinetic, 'rescaled'
      end if

      if (n_Ne > 0) then
         Ndf = 3*n_Ne
         T_kinetic = Ek_Ne*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Ne] = ',T_kinetic
         Ne_vxyz = Ne_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_Ne = kinetic_energy(n_Ne,Ne_vxyz,Ne_mass)
         T_kinetic = Ek_Ne*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Ne] = ',T_kinetic, 'rescaled'
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
         call write_vel(imve,n_He,print_all = .false.)
      end if

      call set_bond_angle_lists()
      call new_tcf(n_He,ntau_msd,nsamp,msd_He,'msd_He.out')
      call new_tcf(n_He,ntau_vcf,nsamp,vcf_He,'vcf_He.out')
      call new_tcf(n_Ne,ntau_msd,nsamp,msd_Ne,'msd_Ne.out')
      call new_tcf(n_Ne,ntau_vcf,nsamp,vcf_Ne,'vcf_Ne.out')

      time = 0.0_wp

      fxyz = 0.0_wp
      He_fxyz = 0.0_wp
      Ne_fxyz = 0.0_wp

      energy = 0.0_wp
!
!     only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_He,He_atom,He_xyz,He_fxyz,ULJEL); print*,'ULJ He  ',sum(ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_Ne,Ne_atom,Ne_xyz,Ne_fxyz,ULJEL); print*,'ULJ Ne  ',sum(ULJEL)
!do i = 1,n_Ne
!   print *,i,ULJEL(i)
!end do
      call force_keating
!     call repul_energy2(1,natom)

!
!======================================================== main loop
      main_loop: do istep = 1, number_moves


      call cpu_time(tp(1))

      coprxyz = rxyz
      call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
      call pbc(rxyz)
      call gas_vverlet_a(dt,He_massi,He_xyz,He_vxyz,He_fxyz,dr_He)
      call pbc(He_xyz)
      r_He_uc = r_He_uc + dr_He

      call gas_vverlet_a(dt,Ne_massi,Ne_xyz,Ne_vxyz,Ne_fxyz,dr_Ne)
      call pbc(Ne_xyz)
      r_Ne_uc = r_Ne_uc + dr_Ne


      ! check if the neighbour list needs to be updated
      do i = 1,natom
         if (CELL(rxyz(i,1:NDIM)) /= CELL(coprxyz(i,1:NDIM))) then
            !! print *,'NEW LIST'
            call NEW_NLIST
            EXIT
         end if
      end do

      fxyz = 0.0_wp
      He_fxyz = 0.0_wp
      Ne_fxyz = 0.0_wp

      energy = 0.0_wp
      call cpu_time(tp(2)); sumtp(1) = sumtp(1) + tp(2)-tp(1)

      call force_keating
!        call repul_energy2(1,natom)
!
      ! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_He,He_atom,He_xyz,He_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_Ne,Ne_atom,Ne_xyz,Ne_fxyz,ULJEL)


      call cpu_time(tp(3)); sumtp(2) = sumtp(2) + tp(3) - tp(2)

!
      call vverlet_b(dt,massi,vxyz,fxyz)
      call gas_vverlet_b(dt,He_massi,He_vxyz,He_fxyz)
      call gas_vverlet_b(dt,Ne_massi,Ne_vxyz,Ne_fxyz)

      time = time + dtreal

      if (istep <= nequib) then
         ! print *,sum(fxyz(1:natom,:),dim=1)
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_He = kinetic_energy(n_He,He_vxyz,He_mass)
         Ek_Ne = kinetic_energy(n_Ne,Ne_vxyz,Ne_mass)

         Ek_tot = Ek + Ek_He + Ek_Ne

         Ndf = 3*(natom + n_He + n_Ne) - 3
         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)

         if (mod(istep,nprint) == 0) then
            write(io_energy,'(8g18.8)') time, EK_tot, energy, & !+ repulenergy*0.5_wp, &
               EK_tot + energy !+ repulenergy*0.5_wp
            write(io_temp,'(8g18.8)') time, T_kinetic
            call flush(io_energy); call flush(io_temp)
         end if
         !if (T_kinetic == 0.0_wp) T_kinetic = 1.0_wp

!        vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
!        He_vxyz = He_vxyz*sqrt(T_kelvin/T_kinetic)
!        Ne_vxyz = Ne_vxyz*sqrt(T_kelvin/T_kinetic)
         call rem_com_momentum(natom,vxyz,mass)
         call rem_com_momentum(n_He,He_vxyz,He_mass)
         call rem_com_momentum(n_Ne,Ne_vxyz,Ne_mass)

         if (n_He > 0) then
            Ndf = 3*n_He
            T_kinetic = Ek_He*2.0_wp/(K_ev*Ndf)
            He_vxyz = He_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_Ne > 0) then
            Ndf = 3*n_Ne
            T_kinetic = Ek_Ne*2.0_wp/(K_ev*Ndf)
            Ne_vxyz = Ne_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (natom > 0) then
         Ndf = 3*natom
         T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
         vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
         end if
      end if

      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then
         call msd_samp(msd_He,r_He_uc)
         call vcf_samp(vcf_He,He_vxyz)
         call msd_samp(msd_Ne,r_Ne_uc)
         call vcf_samp(vcf_Ne,Ne_vxyz)
      end if
      end if


      if (mod(istep,nprint) == 0) then
      if (istep > nequib) then
         ! Kinetic energy of each species
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_He = kinetic_energy(n_He,He_vxyz,He_mass)
         Ek_Ne = kinetic_energy(n_Ne,Ne_vxyz,Ne_mass)
         Ek_tot = Ek + Ek_He+ Ek_Ne
         Ndf = 3*(natom + n_He + n_Ne) - 3
         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)

         write(io_energy,'(8g18.8)') time, Ek_tot, energy, & ! + repulenergy*0.5_wp, &
            Ek_tot + energy !+ repulenergy*0.5_wp
         write(io_temp,'(9g18.8)') time, t_kinetic, &
            2.0_wp*ek/(k_ev*(3*natom)), &
            2.0_wp*ek_He/(k_ev*(3*n_He)), &
            2.0_wp*ek_Ne/(k_ev*(3*n_Ne))
         call flush(io_energy)
         call flush(io_temp)
      end if
      end if

      call cpu_time(tp(4)); sumtp(3) = sumtp(3) + tp(4) - tp(3)
      if (mod(istep,nprint) == 0) then
         write(*,'(i9,8g18.8)') istep,time,tp(2) - tp(1),tp(3) - tp(2),tp(4) - tp(3), tp(4) - tp(1), tp(4) - tp(0)
      end if

      if (make_movie) then
      if (mod(istep,nfreq_movie) == 0) then
         iframe = iframe + 1
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_He,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_He,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
      end if
      end if

! print out the msd
      if (mod( istep-nequib, nsamp*ntau_msd*2 )==0) then
         call  msd_print(msd_He,dtreal)
         call  vcf_print(vcf_He,dtreal)
         call  msd_print(msd_Ne,dtreal)
         call  vcf_print(vcf_Ne,dtreal)
      end if

!===================================================== end main loop
      end do main_loop

      call  msd_print(msd_He,dtreal)
      call  vcf_print(vcf_He,dtreal)
      call  msd_print(msd_Ne,dtreal)
      call  vcf_print(vcf_Ne,dtreal)

      do i = 1,4
         print *,i,sumtp(i)
      end do
      close(io_energy,status='keep')
      close(io_temp,status='keep')

      call get_free_file_unit(io_config_final)
      open(unit=io_config_final,file='final.xyz')
      call write_xyz(io_config_final)
      close(io_config_final,status='KEEP')
END PROGRAM simbox

