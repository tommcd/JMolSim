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
   !
   ! Gets a command argument and assigns it to a variable.
   ! It checks the status of the argument and writes a record to stdout.
   !
   interface get_com_arg
      module procedure get_com_arg_int, get_com_arg_real, &
                       get_com_arg_logical,get_com_arg_string
   end interface
CONTAINS

   SUBROUTINE get_com_arg_int(ia,var,var_name,stat)
      integer,intent(in):: ia
      integer,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE get_com_arg_int

   SUBROUTINE get_com_arg_real(ia,var,var_name,stat)
      integer,intent(in):: ia
      real(wp),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE get_com_arg_real

   SUBROUTINE get_com_arg_logical(ia,var,var_name,stat)
      integer,intent(in):: ia
      logical,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE get_com_arg_logical

   SUBROUTINE get_com_arg_string(ia,var,var_name,stat)
      integer,intent(in):: ia
      character(*),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         write (*,*) 'arg = ', carg(1:nc)
      end if
      var = ''
      var(1:nc) = carg(1:nc)
      write (*,*) var_name,' = ',trim(var)
   END SUBROUTINE get_com_arg_string

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
   public:: rand,get_rand_state,set_rand_state,read_rand_state,write_rand_state
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
   real(wp):: del_rxn, e_activ, etot
   integer:: natom,nseed,natom_max
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
   real(wp),parameter:: bondl_OH = 0.945_wp/angstrom
   real(wp),parameter:: angle_SiOH = (108.5_wp/180.0_wp)*pi
   real(wp),parameter:: bondl_O2  = 0.9699_wp/angstrom
   real(wp),parameter:: bondl_N2  = 1.0464_wp/angstrom
   real(wp),parameter:: bondl_CO2 = 1.161_wp/angstrom
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

!!>include 'constants_mod.f90'

MODULE constants_mod
   USE precision_mod, only: wp
   USE global_vars, only: angstrom,pi
   USE seaton_mod, sigma_OH => sig_O_OH_sil
   implicit none
! TIP3P water
   real(wp),parameter:: bondl_H2O = 0.9572_wp/angstrom
   real(wp),parameter:: angle_H2O = (104.52_wp/180.0_wp)*pi
   integer:: n_CO2
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
   integer,parameter:: iOw = 4
   integer,parameter:: iHw = 5
   integer,parameter:: ntyp = 5
   character(2),parameter:: atom_name(0:ntyp) = (/ 'O ','Si','OH','H ','Ow','Hw' /)
   integer,parameter:: ncmax(0:ntyp) = (/2,4,2,1,1,2/)
   integer,allocatable:: atom(:),CO2_atom(:,:)
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
   real(wp),allocatable:: Ox_xyz(:,:,:),N2_xyz(:,:,:),CO2_xyz(:,:,:)
   real(wp),allocatable:: r_O2_uc(:,:,:),r_N2_uc(:,:,:),r_CO2_uc(:,:,:)
   real(wp),allocatable:: dr_O2(:,:,:),dr_N2(:,:,:),dr_CO2(:,:,:)
   real(wp),allocatable:: rcm_O2_uc(:,:),rcm_N2_uc(:,:),rcm_CO2_uc(:,:)
   real(wp),allocatable:: vcm_O2(:,:),vcm_N2(:,:),vcm_CO2(:,:)
   real(wp),allocatable:: Ox_fxyz(:,:,:),N2_fxyz(:,:,:),CO2_fxyz(:,:,:)
   real(wp),allocatable:: Imige_Ox_xyz(:,:,:),Imige_CO2_xyz(:,:,:),Imige_N2_xyz(:,:,:)
   real(wp),allocatable:: Ox_vxyz(:,:,:),N2_vxyz(:,:,:),CO2_vxyz(:,:,:)
   real(wp),allocatable:: mol_vec(:,:)

   real(wp):: boxl,boxli,boxl2
   interface pbc
      module procedure  pbc_v, pbc_a, pbc_a2
   end interface
CONTAINS

!   PURE SUBROUTINE pbc_v(r)
!      real(wp),intent(inout):: r(:)
!      r(1:3) = r(1:3) - boxl*anint(r(1:3)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_a(r)
!      real(wp),intent(inout):: r(:,:)
!      r(:,1:3) = r(:,1:3) - boxl*anint(r(:,1:3)*boxli)
!   END SUBROUTINE
!
!   PURE SUBROUTINE pbc_a2(r)
!      real(wp),intent(inout):: r(:,:,:)
!      r(:,:,1:3) = r(:,:,1:3) - boxl*anint(r(:,:,1:3)*boxli)
!   END SUBROUTINE
!
   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      if (abs(r(1)) > boxl2) r(1) = r(1) - sign(boxl,r(1))
      if (abs(r(2)) > boxl2) r(2) = r(2) - sign(boxl,r(2))
      if (abs(r(3)) > boxl2) r(3) = r(3) - sign(boxl,r(3))
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      integer:: i
      do i = 1,size(r,1)
         if (abs(r(i,1)) > boxl2) r(i,1) = r(i,1) - sign(boxl,r(i,1))
         if (abs(r(i,2)) > boxl2) r(i,2) = r(i,2) - sign(boxl,r(i,2))
         if (abs(r(i,3)) > boxl2) r(i,3) = r(i,3) - sign(boxl,r(i,3))
      end do
   END SUBROUTINE

   PURE SUBROUTINE pbc_a2(r)
      real(wp),intent(inout):: r(:,:,:)
      integer:: i,j
      do j = 1,size(r,2)
      do i = 1,size(r,1)
         if (abs(r(i,j,1)) > boxl2) r(i,j,1) = r(i,j,1) - sign(boxl,r(i,j,1))
         if (abs(r(i,j,2)) > boxl2) r(i,j,2) = r(i,j,2) - sign(boxl,r(i,j,2))
         if (abs(r(i,j,3)) > boxl2) r(i,j,3) = r(i,j,3) - sign(boxl,r(i,j,3))
      end do
      end do
   END SUBROUTINE


   SUBROUTINE WRITE_XYZ(iu,nattached)
      USE atom_types, only: atom_name,atom
      integer,intent(in):: iu,nattached
      integer:: i
      write(iu,*) natom
      write(iu,*) boxl*Angstrom, nattached
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
   USE seaton_mod, only:bondl_O2,bondl_N2,bondl_CO2
   implicit none
   real(wp),parameter:: KOO = 366.142_wp    ! eV nm^-2
   real(wp),parameter:: KNN = 366.142_wp    ! eV nm^-2
   real(wp),parameter:: KCO2 = 366.142_wp   ! eV nm^-2
   real(wp),parameter:: AOO = bondl_O2     ! nm
   real(wp),parameter:: ANN = bondl_N2     ! nm
   real(wp),parameter:: ACO2 = bondl_CO2    ! nm
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
      USE global_vars, only: natom,angstrom
      USE atom_types, only: atom,atom_name
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 3*nmol + 2
            write(iu,*)'Gases + Silica',boxl
         else
            write(iu,*) 3*nmol + 2
            write(iu,*)'Gases',boxl
         end if
      else
         write(iu,*) 3*nmol + 2
         write(iu,*)'Gases',boxl
      end if
      do i = 1,nmol
         write(iu,110)'CO ',(CO2_xyz(i,1,j)*Angstrom,j = 1,3)
         write(iu,110)'OC ',(CO2_xyz(i,2,j)*Angstrom,j = 1,3)
         write(iu,110)'OC ',(CO2_xyz(i,3,j)*Angstrom,j = 1,3)
      end do
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,'(a,3(1x,f14.6))') atom_name(atom(i)),(rxyz(i,j)*Angstrom,j = 1,3)
         end do
         end if
      end if
      write(iu,110)'Ar ', (/boxl2,boxl2,boxl2/)*Angstrom
      write(iu,110)'Ar ', -(/boxl2,boxl2,boxl2/)*Angstrom
      write(iu,*)
      close(iu,status='keep')
110   format(a,3(1x,f12.6))
   END SUBROUTINE write_box_frame

   SUBROUTINE write_vel(iu,nmol,print_all)
      USE global_vars, only: natom,angstrom
      USE atom_types, only: atom,atom_name
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      integer:: i,j
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 3*nmol
            write(iu,*)'Gases + Silica',boxl
         else
            write(iu,*) 3*nmol
            write(iu,*)'Gases',boxl
         end if
      else
         write(iu,*) 3*nmol + 2
         write(iu,*)'Gases',boxl
      end if
      do i = 1,nmol
         write(iu,*) (CO2_vxyz(i,1,j),j = 1,3)
         write(iu,*) (CO2_vxyz(i,2,j),j = 1,3)
         write(iu,*) (CO2_vxyz(i,3,j),j = 1,3)
      end do
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,*) vxyz(i,1:3)
         end do
         end if
      end if
      write(iu,*)
      close(iu,status='keep')
   END SUBROUTINE write_vel

   SUBROUTINE read_vel(iu,nmol,print_all)
      USE global_vars, only: natom,angstrom
      USE atom_types, only: atom,atom_name
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      real(wp):: rtmp,x,y,z
      integer:: i,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp,rtmp
      do i = 1,nmol
         read(iu,*) x,y,z
         Ox_vxyz(i,1,:) = (/ x,y,z /)
         read(iu,*) x,y,z
         Ox_vxyz(i,2,:) = (/ x,y,z /)
      end do
      do i = 1,nmol
         read(iu,*) x,y,z
         N2_vxyz(i,1,:) = (/ x,y,z /)
         read(iu,*) x,y,z
         N2_vxyz(i,2,:) = (/ x,y,z /)
      end do
      do i = 1,nmol
         read(iu,*) x,y,z
         CO2_vxyz(i,1,:) = (/ x,y,z /)
         read(iu,*) x,y,z
         CO2_vxyz(i,2,:) = (/ x,y,z /)
         read(iu,*) x,y,z
         CO2_vxyz(i,3,:) = (/ x,y,z /)
      end do
      if (present(print_all)) then
      if (print_all) then
         do i = 1,natom
            read(iu,*) x,y,z
            vxyz(i,1:3) = (/ x,y,z /)
         end do
      end if
      end if
      !close(iu,status='keep')
   END SUBROUTINE read_vel

   SUBROUTINE read_box_frame(iu,nmol,print_all)
      USE global_vars, only: natom,angstrom
      USE atom_types, only: atom,atom_name
      USE coordinates_mod
      integer,intent(in):: iu,nmol
      logical,intent(in),optional:: print_all
      character(32):: ctmp
      real(wp):: rtmp,x,y,z
      integer:: i,itmp
!
      read(iu,*) itmp
      read(iu,*) ctmp,rtmp
      do i = 1,nmol
         read(iu,*) ctmp,x,y,z
         Ox_xyz(i,1,:) = (/ x,y,z /)/Angstrom
         read(iu,*) ctmp,x,y,z
         Ox_xyz(i,2,:) = (/ x,y,z /)/Angstrom
      end do
      do i = 1,nmol
         read(iu,*) ctmp,x,y,z
         N2_xyz(i,1,:) = (/ x,y,z /)/Angstrom
         read(iu,*) ctmp,x,y,z
         N2_xyz(i,2,:) = (/ x,y,z /)/Angstrom
      end do
      do i = 1,nmol
         read(iu,*) ctmp,x,y,z
         CO2_xyz(i,1,:) = (/ x,y,z /)/Angstrom
         read(iu,*) ctmp,x,y,z
         CO2_xyz(i,2,:) = (/ x,y,z /)/Angstrom
         read(iu,*) ctmp,x,y,z
         CO2_xyz(i,3,:) = (/ x,y,z /)/Angstrom
      end do
      if (present(print_all)) then
      if (print_all) then
         do i = 1,natom
            read(iu,*) ctmp,x,y,z
            rxyz(i,1:3) = (/ x,y,z /)/Angstrom
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
CONTAINS

   SUBROUTINE COMVEL(natom,n_CO2,n_O2,n_N2,TEMP,vxyz,Ox_vxyz,N2_vxyz,CO2_vxyz)
!   *******************************************************************
!   ** TRANSLATIONAL VELOCITIES FROM MAXWELL-BOLTZMANN DISTRIBUTION  **
!   **                                                               **
!   ** THE DISTRIBUTION IS DETERMINED BY TEMPERATURE AND (unit) MASS.**
!   ** THIS ROUTINE IS GENERAL, AND CAN BE USED FOR ATOMS, LINEAR    **
!   ** MOLECULES, AND NON-LINEAR MOLECULES.                          **
!   **                                                               **
!   ** ROUTINE REFERENCED:                                           **
!   **                                                               **
!   ** real FUNCTION GAUSS ( DUMMY )                                 **
!   **    RETURNS A UNIFORM RANDOM NORMAL VARIATE FROM A             **
!   **    DISTRIBUTION WITH ZERO MEAN AND unit VARIANCE.             **
!   *******************************************************************
!
! INITIAL VELOCITY FOR MOVING PARTICLE
!
      integer,intent(in):: natom,n_CO2,n_O2,n_N2
      real(wp) PI,PI2,RMUL2,IAX1,IAX2,RMO2,IAY1,IAY2,IAZ1,IAZ2
      real(wp):: RIMO1,RNX1,RNX2,RNY1,RNY2,RNZ1,RNZ2
      real(wp):: vxyz(:,:),Ox_vxyz(:,:,:),N2_vxyz(:,:,:),CO2_vxyz(:,:,:),TEMP
      integer:: i,j
      PI = 4.0_wp*ATAN(1.0_wp)
      PI2 = 2.0_wp*PI
      IAX1 = 987321.0_wp
      IAY1 = 199971.0_wp
      IAZ1 = 673923.0_wp
      IAX2 = 786391.0_wp
      IAY2 = 383837.0_wp
      IAZ2 = 921217.0_wp
      RMUL2 = 16807.0_wp
      RMO2 = 2147483647.0_wp
      RIMO1 = 1.0_wp/(2147483648.0_wp)
      do  i = 1,natom + 2*n_O2 + 2*n_N2 + 3*n_CO2
! INITIAL VELOCITY FROM MAXWELLIAN DISTRIBUTION
      IAX1 = MOD(RMUL2*IAX1,RMO2)
      IAX1 = MOD(RMUL2*IAX1,RMO2)
      IAX2 = MOD(RMUL2*IAX2,RMO2)
      IAX2 = MOD(RMUL2*IAX2,RMO2)
      IAY1 = MOD(RMUL2*IAY1,RMO2)
      IAY1 = MOD(RMUL2*IAY1,RMO2)
      IAY2 = MOD(RMUL2*IAY2,RMO2)
      IAY2 = MOD(RMUL2*IAY2,RMO2)
      IAZ1 = MOD(RMUL2*IAZ1,RMO2)
      IAZ1 = MOD(RMUL2*IAZ1,RMO2)
      IAZ2 = MOD(RMUL2*IAZ2,RMO2)
      IAZ2 = MOD(RMUL2*IAZ2,RMO2)
      RNX1 = IAX1*RIMO1
      RNX2 = IAX2*RIMO1
      RNY1 = IAY1*RIMO1
      RNY2 = IAY2*RIMO1
      RNZ1 = IAZ1*RIMO1
      RNZ2 = IAZ2*RIMO1
      vxyz(i,1) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNX1))*COS(PI2*RNX2)
      vxyz(i,2) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNY1))*COS(PI2*RNY2)
      vxyz(i,3) = sqrt(TEMP*K_ev)*SQRT(-2.0_wp*LOG(RNZ1))*COS(PI2*RNZ2)
      end do
         j = 1
      do i = 1,n_O2
       Ox_vxyz(i,1,:) = vxyz(natom + j,:)
       Ox_vxyz(i,2,:) = vxyz(natom + j + 1,:)
        j = j + 2
      end do
         j = 1
      do i = 1,n_N2
       N2_vxyz(i,1,:) = vxyz(natom + 2*n_O2 + j,:)
       N2_vxyz(i,2,:) = vxyz(natom + 2*n_O2 + j + 1,:)
        j = j + 2
      end do
         j = 1
      do i = 1,n_CO2
       CO2_vxyz(i,1,:) = vxyz(natom + 2*n_O2 + 2*n_N2 + j,:)
       CO2_vxyz(i,2,:) = vxyz(natom + 2*n_O2 + 2*n_N2 + j + 1,:)
       CO2_vxyz(i,3,:) = vxyz(natom + 2*n_O2 + 2*n_N2 + j + 2,:)
        j = j + 3
      end do
   END SUBROUTINE

END MODULE

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
      real(wp):: k,dist,rl1(3),r1sq,ktheta,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
          fxyz = 0.0_wp
          energy = 0.0_wp

!     Force routine for Keating bond stretching

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

        energy = energy + 0.5_wp*k*(r1sq - dist)**2

      end do

!     Force routine for Keating bond angle bending

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
        energy = energy + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2

      end do

   END SUBROUTINE
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
      USE coordinates_mod, only: rxyz,boxl2
      implicit none
      PRIVATE
      integer,public:: neigh
      real(wp),private:: delx,dely,delz,delmin
      real(wp),private:: delxi,delyi,delzi
      integer:: ncelx,ncely,ncelz,ncelt
      integer,parameter:: neighmx = 350
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
      CELL = INT((XYZ(1) + boxl2)*delxi) &
           + INT((XYZ(2) + boxl2)*delyi)*ncelx &
           + INT((XYZ(3) + boxl2)*delzi)*ncelx*ncely + 1
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
if (ic < 1) then
 print *,i,rxyz(i,:)
 stop 'NEW_NLIST: ic < 1'
end if
if (ic > ncelt) stop 'NEW_NLIST: ic > ncelt'
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
!               select case(atom(L))
!               case(iOxygenH)
!                  K = proximity(L,1)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  if (proximity(K,3) == M) GOTO 100
!                  if (proximity(K,4) == M) GOTO 100
!                  K = proximity(L,2)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  if (proximity(K,3) == M) GOTO 100
!                  if (proximity(K,4) == M) GOTO 100
!               case(iOxygen)
!                  K = proximity(L,1)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  if (proximity(K,3) == M) GOTO 100
!                  if (proximity(K,4) == M) GOTO 100
!                  K = proximity(L,2)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  if (proximity(K,3) == M) GOTO 100
!                  if (proximity(K,4) == M) GOTO 100
!               case(iSilicon)
!                  K = proximity(L,1)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  K = proximity(L,2)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  K = proximity(L,3)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!                  K = proximity(L,4)
!                  if (K == M) GOTO 100
!                  if (proximity(K,1) == M) GOTO 100
!                  if (proximity(K,2) == M) GOTO 100
!               end select

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

!!>include 'insert_mod.f90'
!!>include 'lj_el_mod_O2_N2_CO2.f90'
!!>include 'force_en_el_mod.f90'
!!>include 'msd_mod.f90'
!!>include 'vcf_mod.f90'

MODULE RDF_MOD
   USE precision_mod, only: wp
   USE global_vars
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

   SUBROUTINE RDF_CALC0
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
            RDF(k) = RDF(k) + 2.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC0

   SUBROUTINE RDF_CALC
      real(wp):: rv(3),r,r2,rmax,zi,zj
      integer:: i,j,k
      DL1 = rc/real(nbin,wp)
      RDF = 0.0_wp
      npairstot = 0.0_wp
      na = 0.0_wp
      do i = nseed + 1,natom
         zi = rxyz(i,3)
         if ( (zi > zu) .or. (zi < zl) ) cycle
         na = na + 1
         do j = nseed + 1,natom
            if (i == j) cycle
            zj = rxyz(j,3)
            if ( (zj > zu + rc) .or. (zj < zl - rc) ) cycle
!            if ( (zi > zu) .or. (zi < zl) ) then
!            if ( (zj > zu) .or. (zj < zl) ) then
!               cycle
!            end if
!            end if
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
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

   SUBROUTINE RDF_CALC_ATOMTYPES
      USE atom_types
      real(wp):: rv(3),r,r2,rmax,rmax2,zi,zj
      integer:: i,j,k,it,ia,ja
      rmax = boxl*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rc/real(nbin,wp)
      RDFT = 0.0_wp
      npairs = 0
      do i = nseed + 1,natom
         ia = atom(i)
         if (ia == iHydrogen) cycle
         zi = rxyz(i,3)
         do j = i + 1,natom
            ja = atom(j)
            if (ja == iHydrogen) cycle
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
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
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

   SUBROUTINE RDF_CALC_ATOMTYPES0
      USE atom_types
      real(wp):: rv(3),r,r2,rmax,rmax2
      integer:: i,j,k,it,ia,ja
      rmax = boxl*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDFT = 0.0_wp
      npairs = 0
      do i = nseed + 1,natom - 1
         ia = atom(i)
         if (ia == iHydrogen) cycle
         do j = i + 1,natom
            ja = atom(j)
            if (ja == iHydrogen) cycle
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
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
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

   SUBROUTINE RDF_PRINT_ATOMTYPES0(io,nattached)
      integer,intent(in):: io,nattached
      integer:: i,ityp
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom - nseed) = ',natom - nseed
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

   SUBROUTINE RDF_PRINT_ATOMTYPES(io,nattached)
      USE atom_types
      integer,intent(in):: io,nattached
      integer:: i,ityp,nati(2),nat
      real(wp):: r,dv,zi,vol
      nat = 0
      nati = 0
      do i = nseed + 1,natom
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
      vol = ((zu - zl)*boxl*boxl)
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom - nseed) = ',natom - nseed
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

   SUBROUTINE RDF_PRINT(io,nattached)
      integer,intent(in):: io,nattached
      integer:: i,nat
      real(wp):: r,dv,zi,vol
      nat = 0
      do i = nseed + 1,natom
         zi = rxyz(i,3)
         if ( (zi <= zu + rc) .and. (zi >= zl - rc) ) then
            nat = nat + 1
         end if
      end do
print * , 'nat = ',nat
print * , 'na  = ',na
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom - nseed) = ',natom - nseed
      write(io,*)'# nat = ',nat
      write(io,*)'# npairstot = ',npairstot
      write(io,*)'# npairstot/(nat - 1) = ',npairstot/(nat - 1)
      rho = nat/((zu - zl)*boxl*boxl)
      vol = ((zu - zl)*boxl*boxl)
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

   SUBROUTINE RDF_PRINT0(io,nattached)
      integer,intent(in):: io,nattached
      integer:: i
      real(wp):: r,dv
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom - nseed) = ',natom - nseed
      rho = (natom - nseed)/(maxval(rxyz(1:natom,3))*boxl*boxl)
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)/(natom - nseed)
      end do
      write(io,'(/)')
      do i = 1,nbin
         r = (i - 0.5_wp)*DL1
         dv = 4.0_wp*pi*(DL1**3)*((i + 1)**3 - i**3)/3.0_wp
         rdf(i) = rdf(i)/(ngr*dv*rho*(natom - nseed))
         write(14,'(f12.6,2x,g16.8)') r,RDF(i)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT0

END MODULE RDF_MOD



PROGRAM SIMBOX
      USE precision_mod
      USE commandline_mod
      USE command_arg_mod
      USE precision_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
!     USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE keating_parameters
      USE nlist_mod
      USE repul_energy_mod
      USE charges_mod
!     USE force_en_el_mod
!     USE insert_mod
!     USE msd_mod
!     USE vcf_mod
      implicit none
      integer,parameter:: number_moves = 1000
      real(wp),parameter:: dtreal = 1.0_wp
      integer,parameter:: n_short = 4
      integer,parameter:: nequib = 200
      integer,parameter:: nsamp = 50
      integer,parameter:: nprint_all = 100
      logical,parameter:: make_movie = .FALSE.
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: mN = 14.00_wp
      real(wp),parameter:: mC = 12.00_wp
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,qj,sumqj,Temp_K
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp)::CL,time,Ek,EkCO2,rad_CO2(3),rr(3),ULJEL,dt_short
      real(wp),allocatable:: massi(:,:),CO2_massi(:,:,:)
      integer:: j,i,m,ii,ic = 0,natomc,istep,stat,nat,nc
      integer:: itmp,i_short,narg
      real(wp),allocatable:: coprxyz(:,:)
      real:: tp(0:10),sumtp(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5
      character(len=132):: infile,carg,ftyp*3
      integer:: io_temp,io_energy,io_config,io_config_final,io_topol
      integer:: io_restart_xyz, io_restart_v,ios,ionum
      logical:: restart
!
      call cpu_time(tp(0))
      open(unit=24,file='final.xyz')
      open(unit=17,file='energy.out')
!
      narg = command_argument_count()
      call get_com_arg(0, carg, "command",stat)
      if (narg < 2) then
         write(*,*)'usage :'
         write(*,*) trim(carg),'  atom_file  T[K]  {xyz_restart  vel_restart} '
         stop
      end if
!
      call get_com_arg(1,infile, "input file",stat)
      io_config = 14
      open(unit=io_config,file=trim(infile))
      call get_com_arg(2,Temp_K, "Temperature [K]",stat)
      call get_com_arg(3, carg, "restart coord file",stat)
      if (stat /= 0) then
         restart = .false.
      else
         restart = .true.
         infile = trim(adjustl(carg))
         write (*,*) 'restarting from ', trim(infile)
         open(unit=34,file=trim(infile))
         call get_com_arg(4, infile, "restart vel file", stat)
         if (stat /= 0) then
            print *,'restart velocity file missing'
            stop
         else
            open(unit=35,file=trim(infile))
         end if
      end if

!
      rad_CO2 = (/sig_C_CO2,sig_O_CO2,sig_O_CO2/)*0.5_wp
      n_CO2 = 0
      T(1) = Temp_K
!
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
      nat = 2
      allocate(CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(r_CO2_uc(n_CO2,3,3))
      allocate(dr_CO2(n_CO2,3,3))
      allocate(rcm_CO2_uc(n_CO2,3))
      allocate(vcm_CO2(n_CO2,3))
      allocate(Imige_CO2_xyz(1:n_CO2,1:3,1:3))

      nc = len_trim(infile)
print *,"'",trim(infile),"'"
print *,nc
      ftyp = infile(nc - 2:nc)
      select case(ftyp)
      case('xyz')
         read (io_config,*) natom
         read (io_config,*) boxl
         io_topol = 15
         open(unit=io_topol,file=trim(infile(1:nc - 3))//'con')
      case('pdb')
         read (io_config,*) c6,natom  ;print *,'natom = ',natom
         read (io_config,*) c6,boxl   ;print *,' boxl = ',boxl
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
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(fxyz(1:natom_max,3),vxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
      allocate(massi(natom_max,3))
      proximity = 0
      allocate(CO2_fxyz(1:n_CO2,1:3,1:3),CO2_vxyz(1:n_CO2,1:3,1:3))
      allocate(CO2_atom(1:n_CO2,3))
      allocate(CO2_gas_charge(n_CO2,3))
      allocate(CO2_massi(1:n_CO2,3,3))
      allocate(mol_vec(3,3))

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
         if (atom(i) == iSilicon) massi(i,:) = 1.0_wp/mSi
         if (atom(i) == iOxygen)  massi(i,:) = 1.0_wp/mOx
         if (atom(i) == iOxygenH)  massi(i,:) = 1.0_wp/mOH
      end do

      select case(ftyp)
      case('xyz')
         read(io_topol,*) itmp
         do i = 1,natom
            read(io_topol,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(io_topol)
      case('pdb')
         do i = 1,natom
            read(io_config,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
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
      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom)
      call NEW_NLIST
!
      if (restart) then
         call read_box_frame(34,nmol = n_CO2,print_all = .true.)
         call read_vel(35,nmol = n_CO2,print_all = .true.)
      end if

      if (make_movie) then ! write out first frame of movie
         iframe = 0
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_CO2,print_all = .true.)
      end if
!========================================================================
!
      call set_bond_angle_lists()
      time = 0.0_wp
      dt = dtreal*dt_fs
      dt_short = dt/n_short

      vxyz = 0.0_wp
      fxyz = 0.0_wp
      energy = 0.0_wp
!
!======================================================== main loop
      main_loop: do istep = 1,number_moves

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi

      fxyz = 0.0_wp      ;call cpu_time(tp(1))
      call force_keating
      call repul_energy2(1,natom)

!=============================== short loop ===
      short_loop: do i_short = 1, n_short
         coprxyz = rxyz
         call vverlet_a(dt_short,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)
!
         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3))) then
               !! print *,'NEW LIST'
               call NEW_NLIST
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         energy = 0.0_wp                           ! ;call cpu_time(tp(4))
         call force_keating                        ! ;call cpu_time(tp(5)); sumtp(4) = sumtp(4) + tp(5) - tp(4)
         call repul_energy2(1,natom)               ! ;call cpu_time(tp(6)); sumtp(5) = sumtp(5) + tp(6) - tp(5)
         call vverlet_b(dt_short,massi,vxyz,fxyz)  ! ;call cpu_time(tp(9)); sumtp(8) = sumtp(8) + tp(9) - tp(8)
!============================end short loop ===
      end do short_loop
      call cpu_time(tp(2)); sumtp(1) = sumtp(1) + tp(2) - tp(1)

      fxyz = 0.0_wp
      vxyz = vxyz + 0.5_wp*dt*fxyz*massi

      time = time + dtreal
      if (istep <= nequib) then
         !
         Ek = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         T(2) = Ek*2.0_wp/( 3.0_wp*K_ev*(natom) )
         !
         if (mod(istep,nsamp) == 0) then
            write(*,*) istep,tp(2) - tp(1),tp(2) - tp(0)
            write(17,'(8g18.8)') time, EK, energy + repulenergy*0.5_wp, Ek + energy + repulenergy*0.5_wp
            call flush(17)
         end if
         vxyz = vxyz*sqrt(T(1)/T(2))
      end if

      if ((mod(istep,nsamp) == 0) .AND. (istep >= nequib)) then
         Ek = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         write (17,'(8g18.8)') time, Ek, Ek + energy + repulenergy*0.5_wp
         call flush(17)
         write(*,*) istep,tp(2) - tp(1),tp(3) - tp(2),tp(3) - tp(0)
         if ( make_movie.and.(mod(istep,nprint_all) == 0) ) then
            iframe = iframe + 1
            call new_frame_file(imve,'frame',iframe)
            call write_box_frame(imve,0,print_all = .TRUE.)
            call new_frame_file(imve,'frame_vel',iframe)
            call write_vel(imve,0,print_all = .TRUE.)
         end if
      end if

!===================================================== end main loop
      end do main_loop

      do i = 1,10
         print *,i,sumtp(i)
      end do
      close(17,status='keep')
      call write_xyz(24,0)
      close(24,status='KEEP')
END PROGRAM simbox

