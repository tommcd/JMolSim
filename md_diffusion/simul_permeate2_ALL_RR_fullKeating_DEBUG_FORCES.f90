
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

END MODULE

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
         allocate(cluster(ncmx,natom_max))
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
   integer:: n_O2,n_N2,n_CO2
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
   integer,allocatable:: atom(:),Ox_atom(:,:),N2_atom(:,:),CO2_atom(:,:)
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
! Parameters for the FULL Keating potential.(KB set)
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
   real(wp),parameter:: KSiSi = 4111.0_wp    ! eV nm^-2
   real(wp),parameter:: KSiO = 68023.0_wp
   real(wp),parameter:: ASiSi = 0.235_wp    ! nm
   real(wp),parameter:: ASiO  = 0.160_wp
   real(wp),parameter:: KSiSiSi =  822.0_wp   ! eV
   real(wp),parameter:: KOSiO  = 25966.0_wp
   real(wp),parameter:: KSiOSi = 11154.0_wp
   real(wp),parameter:: KSiSiO =  4619.96_wp
   real(wp),parameter:: kbond(1:3) = (/ KSiO,KSiSi,KSiO /)
   real(wp),parameter:: abond(1:3) = (/ ASiO,ASiSi,ASiO /)
   real(wp),parameter:: acosang(0:2) = (/ 0.809_wp,1.0_wp/3.0_wp,0.809_wp /)
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
            ctheta(i) = -0.809_wp ! - 1.0_wp
         end if
      end do
!

   END SUBROUTINE set_bond_angle_lists
END MODULE bond_list_mod

!!>include 'box_frames_mod.f90'

MODULE BOX_FRAMES_MOD
   USE precision_mod, only: wp
   USE coordinates_mod
   implicit none
   integer,public:: imve,iframe
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
      integer:: i
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 2*nmol +2*nmol +3*nmol + 2
            write(iu,*)'Gases + Silica',boxl
         else
            write(iu,*) 2*nmol +2*nmol +3*nmol + 2
            write(iu,*)'Gases',boxl
         end if
      else
         write(iu,*) 2*nmol +2*nmol +3*nmol + 2
         write(iu,*)'Gases',boxl
      end if
      do i = 1,nmol
         write(iu,110)'O2 ',Ox_xyz(i,1,:)*Angstrom
         write(iu,110)'O2 ',Ox_xyz(i,2,:)*Angstrom
      end do
      do i = 1,nmol
         write(iu,110)'N2 ',N2_xyz(i,1,:)*Angstrom
         write(iu,110)'N2 ',N2_xyz(i,2,:)*Angstrom
      end do
      do i = 1,nmol
         write(iu,110)'CO ',CO2_xyz(i,1,:)*Angstrom
         write(iu,110)'OC ',CO2_xyz(i,2,:)*Angstrom
         write(iu,110)'OC ',CO2_xyz(i,3,:)*Angstrom
      end do
      if (present(print_all)) then
         if (print_all) then
         do i = 1,natom
            write(iu,'(a,3(1x,f14.6))') atom_name(atom(i)),rxyz(i,1:3)*Angstrom
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
      integer:: i
!
      if (present(print_all)) then
         if (print_all) then
            write(iu,*) natom + 2*nmol +2*nmol +3*nmol
            write(iu,*)'Gases + Silica',boxl
         else
            write(iu,*) 2*nmol +2*nmol +3*nmol
            write(iu,*)'Gases',boxl
         end if
      else
         write(iu,*) 2*nmol +2*nmol +3*nmol + 2
         write(iu,*)'Gases',boxl
      end if
      do i = 1,nmol
         write(iu,*) Ox_vxyz(i,1,:)
         write(iu,*) Ox_vxyz(i,2,:)
      end do
      do i = 1,nmol
         write(iu,*) N2_vxyz(i,1,:)
         write(iu,*) N2_vxyz(i,2,:)
      end do
      do i = 1,nmol
         write(iu,*) CO2_vxyz(i,1,:)
         write(iu,*) CO2_vxyz(i,2,:)
         write(iu,*) CO2_vxyz(i,3,:)
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
      real(wp):: k,dist,rl1(3),r1,ktheta,rl2(3),r2
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
real(wp):: dr1(3),dr2(3),dr3(3),u2,u1,sumdf,ro,dr,df1(3),df2(3),df3(3)
!     fxyz = 0.0_wp
!     energy = 0.0_wp
!     Force routine for Keating bond stretching
      dist = ASiO*ASiO
      do ib = 1,nbondtot
         k = KSiO     !/(10000.0_wp)
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         r1 = dot_product(rl1,rl1)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - 2.0_wp*k*(r1 - dist)*rl1(1:3)
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + 2.0_wp*k*(r1 - dist)*rl1(1:3)
         energy = energy + 0.5_wp*k*(r1 - dist)*(r1 - dist)
      end do

ib = 1
k = KSiO
dist = ASiO**2
ia1 = ibond(1,ib)
ia2 = ibond(2,ib)
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
call pbc(rl1)
r1 = dot_product(rl1,rl1)
df1 = -2.0_wp*k*(r1 - dist)*rl1(1:3)
df2 = + 2.0_wp*k*(r1 - dist)*rl1(1:3)
u1 = 0.5_wp*k*(r1 - dist)*(r1 - dist)
call random_number(dr1)
call random_number(dr2)
dr1 = dr1/1000000.0_wp
dr2 = dr2/1000000.0_wp
write(777,*)
write(777,*),'dr1 = ',dr1
write(777,*),'dr2 = ',dr2
ro = sqrt(r1)
rl1 = rl1 + dr1 - dr2
r1 = dot_product(rl1,rl1)
dr = sqrt(r1) - ro
write(777,*),'dr = ',dr
!df1 =
!df2 =
u2 = 0.5_wp*k*(r1 - dist)*(r1 - dist)
write(777,*),'u1 = ',u1
write(777,*),'u2 = ',u2
write(777,*),'F1.dr1 ',dot_product(df1,dr1)
write(777,*),'F2.dr2 ',dot_product(df2,dr2)
write(777,*),' - F1.dr1 - F2.dr2 = ', -dot_product(df1,dr1) - dot_product(df2,dr2)
write(777,*),'         u2 - u1 = ', u2 - u1
write(777,*),'          F dr = ',sqrt(dot_product(df1,df1))*dr
write(777,*)

      do ib = 1,nbondtot
        k = KSiO
        dist = ASiO
        ia1 = ibond(1,ib)
        ia2 = ibond(2,ib)
        rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
        call pbc(rl1)
        r1 = sqrt(dot_product(rl1,rl1))
        fxyz(ia1,1:3) = fxyz(ia1,1:3) - k*(r1 - dist)*rl1(1:3)/r1
        fxyz(ia2,1:3) = fxyz(ia2,1:3) + k*(r1 - dist)*rl1(1:3)/r1
        energy = energy + 0.5_wp*k*(r1 - dist)**2
      end do

ib = 1
k = KSiO
dist = ASiO
ia1 = ibond(1,ib)
ia2 = ibond(2,ib)
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
call pbc(rl1)
r1 = sqrt(dot_product(rl1,rl1))
df1 = -k*(r1 - dist)*rl1(1:3)/r1
df2 = + k*(r1 - dist)*rl1(1:3)/r1
u1 = 0.5_wp*k*(r1 - dist)**2
call random_number(dr1)
call random_number(dr2)
dr1 = dr1/1000000.0_wp
dr2 = dr2/1000000.0_wp
write(888,*)
write(888,*),'dr1 = ',dr1
write(888,*),'dr2 = ',dr2
ro = r1
rl1 = rl1 + dr1 - dr2
r1 = sqrt(dot_product(rl1,rl1))
dr = r1 - ro
write(888,*),'dr = ',r1 - ro
!df1 = -k*(r1-dist)*rl1(1:3)/r1
!df2 = + k*(r1-dist)*rl1(1:3)/r1
u2 = 0.5_wp*k*(r1 - dist)**2
write(888,*),'u1 = ',u1
write(888,*),'u2 = ',u2
write(888,*),'F1.dr1 ',dot_product(df1,dr1)
write(888,*),'F2.dr2 ',dot_product(df2,dr2)
write(888,*),' - F1.dr1 - F2.dr2 = ', -dot_product(df1,dr1) - dot_product(df2,dr2)
write(888,*),'         u2 - u1 = ', u2 - u1
write(888,*),'          F dr = ',sqrt(dot_product(df1,df1))*dr
write(888,*)

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
         r1 = dot_product(rl1,rl1)
         r2 = dot_product(rl2,rl2)
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1*r2)
         coef = ktheta*(rl1l2 - ctheta(ia)*dist)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - coef*rl2
         fxyz(ia3,1:3) = fxyz(ia3,1:3) - coef*rl1
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + coef*(rl1 + rl2)
         energy = energy + 0.5_wp*ktheta*(rl1l2 - ctheta(ia)*dist)**2
      end do

ia = 1
ktheta = kang(ia)
ia1 = iang(1,ia)
ia2 = iang(2,ia)
ia3 = iang(3,ia)
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
call pbc(rl1)
call pbc(rl2)
r1 = dot_product(rl1,rl1)
r2 = dot_product(rl2,rl2)
rl1l2 = dot_product(rl1,rl2)
cos_theta = (rl1l2)/(r1*r2)
coef = ktheta*(rl1l2 - ctheta(ia)*dist)
df1 = -coef*rl2
df3 = -coef*rl1
df2 = +coef*(rl1 + rl2)
u1 = 0.5_wp*ktheta*(rl1l2 - ctheta(ia)*dist)**2
call random_number(dr1)
call random_number(dr2)
call random_number(dr3)
dr1 = dr1/1000000.0_wp
dr2 = dr2/1000000.0_wp
dr3 = dr3/1000000.0_wp
write(111,*)
write(111,*),'dr1 = ',dr1
write(111,*),'dr2 = ',dr2
write(111,*),'dr3 = ',dr3
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3) + dr1 - dr2
rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3) + dr3 - dr2
call pbc(rl1)
call pbc(rl2)
r1 = dot_product(rl1,rl1)
r2 = dot_product(rl2,rl2)
rl1l2 = dot_product(rl1,rl2)
cos_theta = (rl1l2)/(r1*r2)
coef = ktheta*(rl1l2 - ctheta(ia)*dist)
!df1 = + coef*rl2
!df3 = + coef*rl1
!df2 = -coef*(rl1+rl2)
u2 = 0.5_wp*ktheta*(rl1l2 - ctheta(ia)*dist)**2
write(111,*),'u1 = ',u1
write(111,*),'u2 = ',u2
write(111,*),'F1.dr1 ',dot_product(df1,dr1)
write(111,*),'F2.dr2 ',dot_product(df2,dr2)
write(111,*),'F3.dr3 ',dot_product(df3,dr3)
write(111,*),' - Sum(Fi.dri) = ', -dot_product(df1,dr1) - dot_product(df2,dr2) - dot_product(df3,dr3)
write(111,*),'       u2 - u1 = ', u2 - u1
write(111,*)



!      do ia = 1,nang
!        ktheta = kang(ia)
!        ia1 = iang(1,ia)
!        ia2 = iang(2,ia)
!        ia3 = iang(3,ia)
!        rl1 = rxyz(ia1,1:3)-rxyz(ia2,1:3)
!        rl2 = rxyz(ia3,1:3)-rxyz(ia2,1:3)
!        call pbc(rl1)
!        call pbc(rl2)
!        r1 = sqrt(dot_product(rl1,rl1))
!        r2 = sqrt(dot_product(rl2,rl2))
!        rl1l2 = dot_product(rl1,rl2)
!        cos_theta = (rl1l2)/(r1*r2)
!        coef = ktheta*(cos_theta-ctheta(ia))
!        coef1 = r1**3*r2
!        coef2 = r1*r2
!        coef3 = r1*r2**3
!        fxyz(ia2,1:3) = fxyz(ia2,1:3)-coef*(-(rl1+rl2)/coef2 &
!                                        + rl1l2*rl1/coef1 &
!                                        + rl1l2*rl2/coef3)
!        fxyz(ia1,1:3) = fxyz(ia1,1:3)-coef*(rl2/coef2-rl1l2*rl1/coef1)
!        fxyz(ia3,1:3) = fxyz(ia3,1:3)-coef*(rl1/coef2-rl1l2*rl2/coef3)
!        energy = energy + 0.5_wp*ktheta*(cos_theta-ctheta(ia))**2
!      end do
ia = 1
ktheta = kang(ia)
ia1 = iang(1,ia)
ia2 = iang(2,ia)
ia3 = iang(3,ia)
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
call pbc(rl1)
call pbc(rl2)
r1 = sqrt(dot_product(rl1,rl1))
r2 = sqrt(dot_product(rl2,rl2))
rl1l2 = dot_product(rl1,rl2)
cos_theta = (rl1l2)/(r1*r2)
coef = ktheta*(cos_theta - ctheta(ia))
coef1 = r1**3*r2
coef2 = r1*r2
coef3 = r1*r2**3
df2 = -coef*(-(rl1 + rl2)/coef2  + rl1l2*rl1/coef1 + rl1l2*rl2/coef3)
df1 = -coef*(rl2/coef2 - rl1l2*rl1/coef1)
df3 = -coef*(rl1/coef2 - rl1l2*rl2/coef3)
u1 = 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
call random_number(dr1)
call random_number(dr2)
call random_number(dr3)
dr1 = dr1/1000000.0_wp
dr2 = dr2/1000000.0_wp
dr3 = dr3/1000000.0_wp
write(222,*)
write(222,*),'dr1 = ',dr1
write(222,*),'dr2 = ',dr2
write(222,*),'dr3 = ',dr3
rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3) + dr1 - dr2
rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3) + dr3 - dr2
call pbc(rl1)
call pbc(rl2)
r1 = sqrt(dot_product(rl1,rl1))
r2 = sqrt(dot_product(rl2,rl2))
rl1l2 = dot_product(rl1,rl2)
cos_theta = (rl1l2)/(r1*r2)
coef = ktheta*(cos_theta - ctheta(ia))
coef1 = r1**3*r2
coef2 = r1*r2
coef3 = r1*r2**3
!df2 = -coef*(-(rl1+rl2)/coef2  + rl1l2*rl1/coef1 + rl1l2*rl2/coef3)
!df1 = -coef*(rl2/coef2-rl1l2*rl1/coef1)
!df3 = -coef*(rl1/coef2-rl1l2*rl2/coef3)
u2 = 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
write(222,*),'u1 = ',u1
write(222,*),'u2 = ',u2
write(222,*),'F1.dr1 ',dot_product(df1,dr1)
write(222,*),'F2.dr2 ',dot_product(df2,dr2)
write(222,*),'F3.dr3 ',dot_product(df3,dr3)
write(222,*),' - Sum(Fi.dri) = ', -dot_product(df1,dr1) - dot_product(df2,dr2) - dot_product(df3,dr3)
write(222,*),'       u2 - u1 = ', u2 - u1
write(222,*)

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
      RETURN
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
      RETURN
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
if (ic > ncelt) then
   print *,i,rxyz(i,:)
   stop 'NEW_NLIST: ic > ncelt'
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

   SUBROUTINE insert_mol_N2_O2(n_mol,real_xyz,ntrial,rad1)
      integer,intent(in):: ntrial
      real(wp):: fover = 1.0_wp
      real(wp):: aa(3,3),q(4),x,y,z,real_xyz(:,:,:),rad1(3)
      real(wp):: mol_vec1(3,3)
      integer:: it,i,n_mol,n
!
      mol_loop: do n = 1, n_mol
      insert_loop: do it = 1,ntrial

         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         z = ( 1.0_wp - 2.0_wp*rand() )*boxl2

         !call random_rotation(aa)
         call ran4sph(q)
         call quat2mat(q,aa)
         do i = 1,3
            mol_vec1(:,i) = matmul(aa,mol_vec(:,i))
         end do

         mol_vec1(1,:) = mol_vec1(1,:) + x
         mol_vec1(2,:) = mol_vec1(2,:) + y
         mol_vec1(3,:) = mol_vec1(3,:) + z
         do i = 1,3
            call pbc(mol_vec1(1:3,i))
         end do

! first check if there is a large overlap
         if (.not.overlap(fover,rad1,mol_vec1)) then
            real_xyz(n,1,1) = mol_vec1(1,1)
            real_xyz(n,1,2) = mol_vec1(2,1)
            real_xyz(n,1,3) = mol_vec1(3,1)
            real_xyz(n,2,1) = mol_vec1(1,2)
            real_xyz(n,2,2) = mol_vec1(2,2)
            real_xyz(n,2,3) = mol_vec1(3,2)
            write(*,*) n
            cycle mol_loop
         end if

      end do insert_loop
      end do mol_loop
   END SUBROUTINE insert_mol_N2_O2


   SUBROUTINE insert_mol_CO2(n_mol,real_xyz,ntrial,rad1)
      integer,intent(in):: ntrial
      real(wp):: fover = 1.0_wp
      real(wp):: aa(3,3),q(4),x,y,z,real_xyz(:,:,:),rad1(3)
      real(wp):: mol_vec1(3,3)
      integer:: it,i,n_mol,n

      mol_loop: do n = 1, n_mol
      insert_loop: do it = 1,ntrial
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         z = ( 1.0_wp - 2.0_wp*rand() )*boxl2

         !call random_rotation(aa)
         call ran4sph(q)
         call quat2mat(q,aa)
         do i = 1,3
            mol_vec1(:,i) = matmul(aa,mol_vec(:,i))
         end do

         mol_vec1(1,:) = mol_vec1(1,:) + x
         mol_vec1(2,:) = mol_vec1(2,:) + y
         mol_vec1(3,:) = mol_vec1(3,:) + z
         do i = 1,3
            call pbc(mol_vec1(1:3,i))
         end do

! first check if there is a large overlap
         if (.not.overlap(fover,rad1,mol_vec1)) then
            real_xyz(n,1,1) = mol_vec1(1,1)
            real_xyz(n,1,2) = mol_vec1(2,1)
            real_xyz(n,1,3) = mol_vec1(3,1)
            real_xyz(n,2,1) = mol_vec1(1,2)
            real_xyz(n,2,2) = mol_vec1(2,2)
            real_xyz(n,2,3) = mol_vec1(3,2)
            real_xyz(n,3,1) = mol_vec1(1,3)
            real_xyz(n,3,2) = mol_vec1(2,3)
            real_xyz(n,3,3) = mol_vec1(3,3)
            write(*,*) n
            cycle mol_loop
         end if

      end do insert_loop
      end do mol_loop
   END SUBROUTINE insert_mol_CO2


   FUNCTION overlap(fr,rad,mol_vec1)
      logical:: overlap
      integer:: k,jj,j,ic,nc
      real(wp):: dr(3),rad(3),fr,mol_vec1(3,3)
      overlap = .false.
      probe_atom_loop: do k = 1,3
         ic = CELL( mol_vec1(1:3,k) )
         call NEIGCELL(ic,2,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0) cycle cell_loop
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               if (atom(j) /= iSilicon) then
               dr(1:3) = mol_vec1(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               if (dot_product(dr,dr) < ( fr*(rad(k) + sigi2(atom(j)) ) )**2) then
                  overlap = .true.
                  RETURN
               end if
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
      end do probe_atom_loop
   END FUNCTION


   FUNCTION overlap2(fr,rad,mol_vec1)
      logical:: overlap2
      integer:: k,j
      real(wp):: dr(3),rad(3),fr,mol_vec1(3,3)
      overlap2 = .false.
      probe_atom_loop: do k = 1,3
         do j = 1,natom
            if (atom(j) == iSilicon) cycle
            dr(1:3) = mol_vec1(1:3,k) - rxyz(j,1:3)
            call pbc(dr)
            if (dot_product(dr,dr) < ( fr*(rad(k) + sigi2(atom(j)) ) )**2) then
               overlap2 = .true.
               RETURN
            end if
         end do
      end do probe_atom_loop
   END FUNCTION

END MODULE insert_mod

!!>include 'lj_el_mod_O2_N2_CO2.f90'

MODULE Lj_el_mod
    USE precision_mod
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
!   real(wp),parameter:: rcut_1 = 3.32_wp/angstrom
!   real(wp),parameter:: rcut_1_2 = rcut_1**2
!   real(wp),parameter:: rcut_1_6 = 1.0_wp/(rcut_1_2*rcut_1_2*rcut_1_2)
!   real(wp),parameter:: rcut_1_7 = 1.0_wp/(((rcut_1**3)**2)*rcut_1)
!   real(wp),parameter:: rcut_1_13 = rcut_1_7*rcut_1_7*rcut_1
    real(wp),parameter:: theta_CO2 = 3.1415926535
    real(wp),parameter:: K_CO2 = 5.331

CONTAINS

   SUBROUTINE LJ_INIT
      real(wp):: eij,rstij
      integer:: i,j
      allocate( aij(0:ntyplj,0:ntyplj) )
      allocate( bij(0:ntyplj,0:ntyplj) )
      print *,'ntyplj = ',ntyplj
      allocate( vdwcut(0:ntyplj,0:ntyplj) )
      allocate( vdwsf(0:ntyplj,0:ntyplj) )
      allocate( vdwcut_mod(0:ntyplj,0:ntyplj) )
      allocate( vdwsf_mod(0:ntyplj,0:ntyplj) )
!
      do j = 0,ntyplj    ! lj aij,bij parameters & shifted terms
      do i = 0,ntyplj
         eij = sqrt(epsi(i)*epsi(j))
         rstij = (sigi(i) + sigi(j))*0.5
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


   SUBROUTINE FORCE_ENERGY_LJ_EL_O2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(3),df(3),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         Imige_Ox_xyz(M,1,:) = Ox_xyz(M,1,:)
         rr(1:3) = Ox_xyz(M,2,:) - Ox_xyz(M,1,:)
         call pbc(rr)
         Imige_Ox_xyz(M,2,:) = Ox_xyz(M,1,:) + rr(:)
         Imige_Ox_xyz(M,3,:) = (Imige_Ox_xyz(M,1,:) + Imige_Ox_xyz(M,2,:))/2.0_wp

       do n = 1,2
         rr(1:3) = Ox_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2*r2*r2)
            rr13 = rr6*rr6/r1
            rr7 = rr6/r1
            dele = rr6*(aij(atom(L),Ox_atom(M,n))*rr6 - bij(atom(L),Ox_atom(M,n))) &
                 - vdwcut(atom(L),Ox_atom(M,n)) - vdwsf(atom(L),Ox_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),Ox_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),Ox_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            Ox_fxyz(M,n,1:3) = Ox_fxyz(M,n,1:3) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,Ox_atom(M,n))*rr6-bij(0,Ox_atom(M,n))) &
!                 - vdwcut_mod(0,Ox_atom(M,n)) - vdwsf_mod(0,Ox_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,Ox_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,Ox_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:3)/r1
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            Ox_fxyz(M,n,1:3) = Ox_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + dele
!         end if

!         rr(1:3) = Imige_Ox_xyz(M,n,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = rr(1)**2+rr(2)**2+rr(3)**2
!         r1 = sqrt(r2)
!            df = -charge(L)*Ox_gas_charge(M,n)*rr(1:3)/(r2*r1)
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            Ox_fxyz(M,n,1:3) = Ox_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,n)/r1
!!

      end do
!
!         rr(1:3) = Imige_Ox_xyz(M,3,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = rr(1)**2+rr(2)**2+rr(3)**2
!         r1 = sqrt(r2)
!         df = -charge(L)*Ox_gas_charge(M,3)*rr(1:3)/(r2*r1)
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            Ox_fxyz(M,1,1:3) = Ox_fxyz(M,1,1:3) - 0.5_wp*df
!            Ox_fxyz(M,2,1:3) = Ox_fxyz(M,2,1:3) - 0.5_wp*df
!         ULJEL = ULJEL + charge(L)*Ox_gas_charge(M,3)/r1
!
!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,Ox_atom(M,2))*rr6-bij(0,Ox_atom(M,2))) &
!                 - vdwcut_mod(0,Ox_atom(M,2)) - vdwsf_mod(0,Ox_atom(M,2))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,Ox_atom(M,2))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,Ox_atom(M,2))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:3)/r1
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            Ox_fxyz(M,1,1:3) = Ox_fxyz(M,1,1:3) - 0.5_wp*df
!            Ox_fxyz(M,2,1:3) = Ox_fxyz(M,2,1:3) - 0.5_wp*df
!            ULJEL = ULJEL + dele
!         end if

      end do atom_loop
      end do outer_atom_loop
      energy = energy + ULJEL
   END SUBROUTINE


   SUBROUTINE FORCE_ENERGY_LJ_EL_N2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(3),df(3),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         Imige_N2_xyz(M,1,:) = N2_xyz(M,1,:)
         rr(1:3) = N2_xyz(M,2,:) - N2_xyz(M,1,:)
         call pbc(rr)
         Imige_N2_xyz(M,2,:) = N2_xyz(M,1,:) + rr(:)
         Imige_N2_xyz(M,3,:) = (Imige_N2_xyz(M,1,:) + Imige_N2_xyz(M,2,:))/2.0_wp

       do n = 1,2
         rr(1:3) = N2_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            dele = rr6*(aij(atom(L),N2_atom(M,n))*rr6 - bij(atom(L),N2_atom(M,n))) &
                 - vdwcut(atom(L),N2_atom(M,n)) - vdwsf(atom(L),N2_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),N2_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),N2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            N2_fxyz(M,n,1:3) = N2_fxyz(M,n,1:3) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,N2_atom(M,n))*rr6-bij(0,N2_atom(M,n))) &
!                 - vdwcut_mod(0,N2_atom(M,n)) - vdwsf_mod(0,N2_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,N2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,N2_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:3)/r1
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            N2_fxyz(M,n,1:3) = N2_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + dele
!         end if
!
!         rr(1:3) = Imige_N2_xyz(M,n,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = rr(1)**2+rr(2)**2+rr(3)**2
!         r1 = sqrt(r2)
!            df = -charge(L)*N2_gas_charge(M,n)*rr(1:3)/(r2*r1)
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            N2_fxyz(M,n,1:3) = N2_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + charge(L)*N2_gas_charge(M,n)/r1
      end do

!         rr(1:3) = Imige_N2_xyz(M,3,:) - rxyz(L,:)
!         call pbc(rr)
!         r2 = rr(1)**2+rr(2)**2+rr(3)**2
!         r1 = sqrt(r2)
!         df = -charge(L)*N2_gas_charge(M,3)*rr(1:3)/(r2*r1)
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            N2_fxyz(M,1,1:3) = N2_fxyz(M,1,1:3) - 0.5_wp*df
!            N2_fxyz(M,2,1:3) = N2_fxyz(M,2,1:3) - 0.5_wp*df
!         ULJEL = ULJEL + charge(L)*N2_gas_charge(M,3)/r1
!
!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,N2_atom(M,2))*rr6-bij(0,N2_atom(M,2))) &
!                 - vdwcut_mod(0,N2_atom(M,2)) - vdwsf_mod(0,N2_atom(M,2))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,N2_atom(M,2))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,N2_atom(M,2))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:3)/r1
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            N2_fxyz(M,1,1:3) = N2_fxyz(M,1,1:3) - 0.5_wp*df
!            N2_fxyz(M,2,1:3) = N2_fxyz(M,2,1:3) - 0.5_wp*df
!            ULJEL = ULJEL + dele
!         end if

      end do atom_loop
      end do outer_atom_loop
      energy = energy + ULJEL
   END SUBROUTINE

    SUBROUTINE FORCE_ENERGY_LJ_EL_CO2(atomfirst,atomlast,sysfirst,syslast,ULJEL)
      integer,intent(in):: atomfirst,atomlast,sysfirst,syslast
      real(wp),intent(out):: ULJEL
      real(wp):: rr(3),df(3),r2,r1,rr6,dele,rr7,rr13,rcut_1_7,rcut_1_13
      integer:: L,M,n
      ULJEL = 0.0_wp
      outer_atom_loop: do L = atomfirst,atomlast
      atom_loop: do M = sysfirst,syslast
         Imige_CO2_xyz(M,1,:) = CO2_xyz(M,1,:)
         rr(1:3) = CO2_xyz(M,2,:) - CO2_xyz(M,1,:)
         call pbc(rr)
         Imige_CO2_xyz(M,2,:) = CO2_xyz(M,1,:) + rr(:)
         rr(1:3) = CO2_xyz(M,3,:) - CO2_xyz(M,1,:)
         call pbc(rr)
         Imige_CO2_xyz(M,3,:) = CO2_xyz(M,1,:) + rr(:)

         do n = 1,3
         rr(1:3) = CO2_xyz(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         r1 = sqrt(r2)

         if (r2 < rcut2) then
            rr6 = 1.0_wp/(r2**3)
            rr13 = rr6*rr6/r1
            rr7 = 1.0_wp/(r1*(r2**3))
            dele = rr6*(aij(atom(L),CO2_atom(M,n))*rr6 - bij(atom(L),CO2_atom(M,n))) &
                 - vdwcut(atom(L),CO2_atom(M,n)) - vdwsf(atom(L),CO2_atom(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),CO2_atom(M,n))*(1.0_wp*rr7 - 1.0_wp*rcut7)&
                      + 12.0_wp*aij(atom(L),CO2_atom(M,n))*(1.0_wp*rcut13 - 1.0_wp*rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            CO2_fxyz(M,n,1:3) = CO2_fxyz(M,n,1:3) - df
            ULJEL = ULJEL + dele
         end if

!         if ((r2 < rcut_1_2).AND.(atom(L) == 1)) then
!            rr6 = 1.0_wp/(r2**3)
!            rr13 = rr6*rr6/r1
!            rr7 = 1.0_wp/(r1*(r2**3))
!            dele = rr6*(aij(0,CO2_atom(M,n))*rr6-bij(0,CO2_atom(M,n))) &
!                 - vdwcut_mod(0,CO2_atom(M,n)) - vdwsf_mod(0,CO2_atom(M,n))*(r1-rcut_1)
!            df = (6.0_wp*bij(0,CO2_atom(M,n))*(1.0_wp*rr7-1.0_wp*rcut_1_7)&
!                      + 12.0_wp*aij(0,CO2_atom(M,n))*(1.0_wp*rcut_1_13 - 1.0_wp*rr13))*rr(1:3)/r1
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            CO2_fxyz(M,n,1:3) = CO2_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + dele
!         end if
!
!            rr(1:3) = Imige_CO2_xyz(M,n,:) - rxyz(L,:)
!            call pbc(rr)
!            r2 = rr(1)**2+rr(2)**2+rr(3)**2
!            r1 = sqrt(r2)
!
!            df = -charge(L)*CO2_gas_charge(M,n)*rr(1:3)/(r2*r1)
!            fxyz(L,1:3) = fxyz(L,1:3) + df
!            CO2_fxyz(M,n,1:3) = CO2_fxyz(M,n,1:3) - df
!            ULJEL = ULJEL + charge(L)*CO2_gas_charge(M,n)/r1


         end do

      end do atom_loop
      end do outer_atom_loop
      energy = energy + ULJEL
   END SUBROUTINE


   SUBROUTINE O2_stretching
!     Force routine for Keating bond stretching
      real(wp):: rl1(3),df(3),r1sq
      integer:: i
      do i = 1,n_O2
        rl1 = Ox_xyz(i,1,:) - Ox_xyz(i,2,:)
        call pbc(rl1)
        r1sq = sqrt(dot_product(rl1,rl1))
        Ox_fxyz(i,1,1:3) = Ox_fxyz(i,1,1:3) - KOO*(r1sq - AOO)*rl1(1:3)/r1sq
        Ox_fxyz(i,2,1:3) = Ox_fxyz(i,2,1:3) + KOO*(r1sq - AOO)*rl1(1:3)/r1sq
        energy = energy + 0.5_wp*KOO*(r1sq - AOO)**2
      end do
    END SUBROUTINE O2_stretching

   SUBROUTINE N2_stretching
!     Force routine for Keating bond stretching
      real(wp):: r1(3),df(3),r1sq
      integer:: i
      do i = 1,n_N2
        r1 = N2_xyz(i,1,:) - N2_xyz(i,2,:)
        call pbc(r1)
        r1sq = sqrt(dot_product(r1,r1))
        N2_fxyz(i,1,1:3) = N2_fxyz(i,1,1:3) - KNN*(r1sq - ANN)*r1(1:3)/r1sq
        N2_fxyz(i,2,1:3) = N2_fxyz(i,2,1:3) + KNN*(r1sq - ANN)*r1(1:3)/r1sq
        energy = energy + 0.5_wp*KNN*(r1sq - ANN)**2
      end do
    END SUBROUTINE


   SUBROUTINE CO2_stretching
!     Force routine for Keating bond stretching
      real(wp):: rl1(3),rl2(3),df1(3),df2(3),r1sq,r2sq,dist
      integer:: i
      dist = ACO2
      do i = 1,n_CO2
        rl1 = CO2_xyz(i,2,:) - CO2_xyz(i,1,:)
        call pbc(rl1)
        r1sq = sqrt(dot_product(rl1,rl1))
        df1 = KCO2*(r1sq - dist)*rl1(1:3)/r1sq
        CO2_fxyz(i,2,1:3) = CO2_fxyz(i,2,1:3) - df1
        CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3) + df1
        energy = energy + 0.5_wp*KCO2*(r1sq - dist)**2
!
        rl2 = CO2_xyz(i,3,:) - CO2_xyz(i,1,:)
        call pbc(rl2)
        r2sq = sqrt(dot_product(rl2,rl2))
        df2 = KCO2*(r2sq - dist)*rl2(1:3)/r2sq
        CO2_fxyz(i,3,1:3) = CO2_fxyz(i,3,1:3) - df2
        CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3) + df2
        energy = energy + 0.5_wp*KCO2*(r2sq - dist)**2
      end do
    END SUBROUTINE


    SUBROUTINE CO2_angle_bending
!     Force routine for CO2 bond angle bending
      real(wp):: rl1(3),rl2(3),df1(3),df2(3),df3(3)
      real(wp):: r1sq,r2sq,rl1l2,cos_theta,theta,coef,coef1,coef2,coef3
      integer:: i
      do i = 1,n_CO2
         rl1 = CO2_xyz(i,2,1:3) - CO2_xyz(i,1,1:3)
         rl2 = CO2_xyz(i,3,1:3) - CO2_xyz(i,1,1:3)
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
!         CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3)-coef*(-(rl1+rl2)/coef2 &
!                           + rl1l2*rl1/coef1 &
!                           + rl1l2*rl2/coef3)
!         CO2_fxyz(i,2,1:3) = CO2_fxyz(i,2,1:3)-coef*(rl2/coef2-rl1l2*rl1/coef1)
!         CO2_fxyz(i,3,1:3) = CO2_fxyz(i,3,1:3)-coef*(rl1/coef2-rl1l2*rl2/coef3)
         df1(1:3) =  coef*( (rl1 + rl2)/coef2 - rl1l2*rl1/coef1 - rl1l2*rl2/coef3 )
         df2(1:3) =  coef*( rl1l2*rl1/coef1 - rl2/coef2 )
         df3(1:3) =  coef*( rl1l2*rl2/coef3 - rl1/coef2 )

         CO2_fxyz(i,1,1:3) = CO2_fxyz(i,1,1:3) + df1
         CO2_fxyz(i,2,1:3) = CO2_fxyz(i,2,1:3) + df2
         CO2_fxyz(i,3,1:3) = CO2_fxyz(i,3,1:3) + df3
         end if
         energy = energy + 0.5_wp*K_CO2*(theta - theta_CO2)**2
      end do
   END SUBROUTINE


   SUBROUTINE FORCE_ENERGY_LJ(sysfirst,syslast,gasfirst,gaslast,nat_gas, &
                             atom_gas,rxyz_gas,fxyz_gas,ULJEL)
      integer,intent(in):: sysfirst,syslast,gasfirst,gaslast,nat_gas
      integer,intent(in):: atom_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(3),df(3),r2,r1,rr6,dele,rr7,rr13
      integer:: L,M,n
      ULJEL = 0.0_wp
      solid_loop: do L = sysfirst,syslast
      mol_loop: do M = gasfirst,gaslast
      atom_loop: do n = 1,nat_gas
         rr(1:3) = rxyz_gas(M,n,:) - rxyz(L,:)
         call pbc(rr)
         r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
         if (r2 < rcut2) then
            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7
            dele = rr6*(aij(atom(L),atom_gas(M,n))*rr6 - bij(atom(L),atom_gas(M,n))) &
                 -   vdwcut(atom(L),atom_gas(M,n)) - vdwsf(atom(L),atom_gas(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),atom_gas(M,n))*(rr7 - rcut7) &
               + 12.0_wp*aij(atom(L),atom_gas(M,n))*(rcut13 - rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            fxyz_gas(M,n,1:3) = fxyz_gas(M,n,1:3) - df
            ULJEL = ULJEL + dele
         end if
      end do atom_loop
      end do mol_loop
      end do solid_loop
      energy = energy + ULJEL
   END SUBROUTINE

   SUBROUTINE FORCE_ENERGY_LJ_NLIST(gasfirst,gaslast,nat_gas, &
              atom_gas,rxyz_gas,fxyz_gas,ULJEL)

      integer,intent(in):: gasfirst,gaslast,nat_gas
      integer,intent(in):: atom_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: ULJEL
      real(wp):: rr(3),df(3),r2,r1,rr6,dele,rr7,rr13
      integer,parameter:: nlm = 5
      integer:: L,M,n,ic,jj,nc,ncell((2*nlm + 1)**3),neigh,nl
      nl = nlayers_ll(rcut - tiny(1.0))
      if (nl > nlm) STOP 'nl > nlm '
      ULJEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
      atom_loop: do n = 1,nat_gas
         ic = CELL( rxyz_gas(M,n,:) )
         call NEIGCELL(ic,nl,neigh,ncell)
         cell_loop: do jj = 1,neigh
         nc = ncell(jj)
         if (nc == 0) cycle cell_loop
         L = HOC(nc)
         cell_atom_loop: do while (L /= 0)
            rr(1:3) = rxyz_gas(M,n,:) - rxyz(L,:)
            call pbc(rr)
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            if (r2 < rcut2) then
            r1 = sqrt(r2)
            rr6 = 1.0_wp/(r2**3)
            rr7 = rr6/r1
            rr13 = rr6*rr7
            dele = rr6*(aij(atom(L),atom_gas(M,n))*rr6 - bij(atom(L),atom_gas(M,n))) &
                 -   vdwcut(atom(L),atom_gas(M,n)) - vdwsf(atom(L),atom_gas(M,n))*(r1 - rcut)
            df = (6.0_wp*bij(atom(L),atom_gas(M,n))*(rr7 - rcut7) &
               + 12.0_wp*aij(atom(L),atom_gas(M,n))*(rcut13 - rr13))*rr(1:3)/r1
            fxyz(L,1:3) = fxyz(L,1:3) + df
            fxyz_gas(M,n,1:3) = fxyz_gas(M,n,1:3) - df
            ULJEL = ULJEL + dele
            end if
         L = LL(L)
         end do cell_atom_loop
         end do cell_loop
      end do atom_loop
      end do mol_loop
      energy = energy + ULJEL
   END SUBROUTINE

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
      real(wp):: rgas(3),rr(3),df(3),r2,r1,du
      integer:: L,M,n
      UEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
         atom_loop: do n = 1,nat_gas
            rr = rxyz_gas(M,n,:) - rxyz_gas(M,1,:)
            call pbc(rr)
            rgas = rxyz_gas(M,1,:) + rr
            solid_loop: do L = sysfirst,syslast
               rr(1:3) = rgas - rxyz(L,:)
               call pbc(rr)
               r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
               r1 = sqrt(r2)
               du = charge(L)*charge_gas(M,n)/r1
               df = du*rr/r2
               fxyz(L,1:3) = fxyz(L,1:3) - df
               fxyz_gas(M,n,1:3) = fxyz_gas(M,n,1:3) + df
               Uel = Uel + du
            end do solid_loop
         end do atom_loop
      end do mol_loop
      energy = energy + UEL
   END SUBROUTINE


   SUBROUTINE FORCE_ENERGY_EL_X2_dummy(sysfirst,syslast,gasfirst,gaslast, &
                                       charge_gas,rxyz_gas,fxyz_gas,UEL)
      integer,intent(in):: sysfirst,syslast,gasfirst,gaslast
      real(wp),intent(in):: charge_gas(:,:)
      real(wp),intent(in):: rxyz_gas(:,:,:)
      real(wp),intent(inout):: fxyz_gas(:,:,:)
      real(wp),intent(out):: UEL
      real(wp):: rgas(3,3),rr(3),df(3),r2,r1,du
      integer:: L,M,n
      UEL = 0.0_wp
      mol_loop: do M = gasfirst,gaslast
         rgas(1,:) = rxyz_gas(M,1,:)
         rr(1:3) = rxyz_gas(M,2,:) - rxyz_gas(M,1,:)
         call pbc(rr)
         rgas(2,:) = rxyz_gas(M,1,:) + rr(:)
         rgas(3,:) = (rgas(1,:) + rgas(2,:))*0.5_wp
         solid_loop: do L = sysfirst,syslast
            do n = 1,2
               rr(1:3) = rgas(n,:) - rxyz(L,:)
               call pbc(rr)
               r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
               r1 = sqrt(r2)
               du = charge(L)*charge_gas(M,n)/r1
               df = du*rr/r2
               fxyz(L,1:3) = fxyz(L,1:3) - df
               fxyz_gas(M,n,1:3) = fxyz_gas(M,n,1:3) + df
               Uel = Uel + du
            end do
            rr(1:3) = rgas(3,:) - rxyz(L,:)
            call pbc(rr)
            r2 = rr(1)**2 + rr(2)**2 + rr(3)**2
            r1 = sqrt(r2)
            du = charge(L)*charge_gas(M,3)/r1
            df = du*rr/r2
            fxyz(L,1:3) = fxyz(L,1:3) - df
            fxyz_gas(M,1,1:3) = fxyz_gas(M,1,1:3) + 0.5_wp*df
            fxyz_gas(M,2,1:3) = fxyz_gas(M,2,1:3) + 0.5_wp*df
            UEL = Uel + du
         end do solid_loop
      end do mol_loop
      energy = energy + UEL
   END SUBROUTINE FORCE_ENERGY_EL_X2_dummy

END MODULE

!!>include 'msd_mod.f90'

MODULE msd_mod
   USE precision_mod
   implicit none
   private
   save
   integer:: icount = 0,isamp = 0,ntau_msd
   real(wp),allocatable:: rOx(:,:,:),rN2(:,:,:),rCO2(:,:,:)
   real(wp),allocatable:: rac1Ox(:,:),rac1N2(:,:),rac1CO2(:,:)
   real(wp),allocatable:: rac1Ox_av(:),rac1N2_av(:),rac1CO2_av(:)
   public:: msd_init,msd_samp,msd_print,ntau_msd
!
CONTAINS

   SUBROUTINE msd_init(nmmol)
      integer,intent(in):: nmmol
      allocate( rOx(nmmol,ntau_msd,3),rN2(nmmol,ntau_msd,3),rCO2(nmmol,ntau_msd,3) )
      allocate( rac1Ox(nmmol,ntau_msd),rac1N2(nmmol,ntau_msd),rac1CO2(nmmol,ntau_msd) )
      allocate( rac1Ox_av(ntau_msd),rac1N2_av(ntau_msd),rac1CO2_av(ntau_msd) )
      rOx(:,:,:) = 0.0
      rN2(:,:,:) = 0.0
      rCO2(:,:,:) = 0.0
      rac1Ox(:,:) = 0.0
      rac1N2(:,:) = 0.0
      rac1CO2(:,:) = 0.0
      RETURN
   END SUBROUTINE msd_init

   SUBROUTINE msd_samp(nmmol,rcmO2,rcmN2,rcmCO2)
      integer,intent(in):: nmmol
      real(wp),intent(in):: rcmO2(:,:),rcmN2(:,:),rcmCO2(:,:)
      integer:: i,l,m
      integer,save:: itau = 0
!
      icount = icount + 1
      itau = itau + 1
      if (itau > ntau_msd) itau = itau - ntau_msd

      forall(i = 1:nmmol)
         rOx(i,itau,:) = rcmO2(i,:)
         rN2(i,itau,:) = rcmN2(i,:)
         rCO2(i,itau,:) = rcmCO2(i,:)
      end forall

      if (icount <= ntau_msd) RETURN
      isamp = isamp + 1
      do l = 1,ntau_msd
         m = itau - l + 1
         if (m <= 0) m = m + ntau_msd
         forall(i = 1:nmmol)
            rac1Ox(i,l) = rac1Ox(i,l) + ( rOx(i,itau,1) - rOx(i,m,1) )**2 &
                                      + ( rOx(i,itau,2) - rOx(i,m,2) )**2 &
                                      + ( rOx(i,itau,3) - rOx(i,m,3) )**2
            rac1N2(i,l) = rac1N2(i,l) + ( rN2(i,itau,1) - rN2(i,m,1) )**2 &
                                      + ( rN2(i,itau,2) - rN2(i,m,2) )**2 &
                                      + ( rN2(i,itau,3) - rN2(i,m,3) )**2
            rac1CO2(i,l) = rac1CO2(i,l) + ( rCO2(i,itau,1) - rCO2(i,m,1) )**2 &
                                        + ( rCO2(i,itau,2) - rCO2(i,m,2) )**2 &
                                        + ( rCO2(i,itau,3) - rCO2(i,m,3) )**2
         end forall
      end do
      RETURN
   END SUBROUTINE msd_samp


   SUBROUTINE msd_print(nmmol,nfreq,dtfs)
      integer,intent(in):: nmmol,nfreq
      real(wp):: dtfs
      real(wp):: rnmal,tau
      integer:: j
      if (isamp <= 0) RETURN
      rnmal = 1.0_wp/isamp
! sum individual molecule tracer MSDs
      rac1Ox_av(:) = 0.0
      rac1N2_av(:) = 0.0
      rac1CO2_av(:) = 0.0
      do j = 1,nmmol
         rac1Ox_av(:) = rac1Ox_av(:) + rac1Ox(j,:)
         rac1N2_av(:) = rac1N2_av(:) + rac1N2(j,:)
         rac1CO2_av(:) = rac1CO2_av(:) + rac1CO2(j,:)
      end do
      rac1Ox_av(:) = rac1Ox_av(:)*rnmal/nmmol
      rac1N2_av(:) = rac1N2_av(:)*rnmal/nmmol
      rac1CO2_av(:) = rac1CO2_av(:)*rnmal/nmmol

      do j = 1,ntau_msd
         tau = (j - 1)*nfreq*dtfs*0.001_wp
         write(1001,*) tau,rac1Ox_av(j)
         write(1002,*) tau,rac1N2_av(j)
         write(1003,*) tau,rac1CO2_av(j)
      end do
      write(1001,'(/)')
      write(1002,'(/)')
      write(1003,'(/)')
   END SUBROUTINE msd_print

END MODULE msd_mod

!!>include 'vcf_mod.f90'

MODULE vcf_mod
   USE precision_mod, only: wp
   implicit none
   private
   save
   integer:: icount = 0,isamp = 0,ntau_vcf
   real(wp),allocatable:: rOx(:,:,:),rN2(:,:,:),rCO2(:,:,:)
   real(wp),allocatable:: rac1Ox(:,:),rac1N2(:,:),rac1CO2(:,:)
   real(wp),allocatable:: rac1Ox_av(:),rac1N2_av(:),rac1CO2_av(:)
   public:: vcf_init,vcf_samp,vcf_print,ntau_vcf
!
CONTAINS

   SUBROUTINE vcf_init(nmmol)
      USE files_mod
      implicit none
      integer,intent(in):: nmmol
      allocate( rOx(nmmol,ntau_vcf,3),rN2(nmmol,ntau_vcf,3),rCO2(nmmol,ntau_vcf,3) )
      allocate( rac1Ox(nmmol,ntau_vcf),rac1N2(nmmol,ntau_vcf),rac1CO2(nmmol,ntau_vcf) )
      allocate( rac1Ox_av(ntau_vcf),rac1N2_av(ntau_vcf),rac1CO2_av(ntau_vcf) )
      rOx(:,:,:) = 0.0
      rN2(:,:,:) = 0.0
      rCO2(:,:,:) = 0.0
      rac1Ox(:,:) = 0.0
      rac1N2(:,:) = 0.0
      rac1CO2(:,:) = 0.0
      RETURN
   END SUBROUTINE vcf_init

   SUBROUTINE vcf_samp(nmmol,rcmO2,rcmN2,rcmCO2)
      implicit none
      integer,intent(in):: nmmol
      real(wp),intent(in):: rcmO2(:,:),rcmN2(:,:),rcmCO2(:,:)
      integer:: i,j,l,m
      integer,save:: itau = 0

      icount = icount + 1
      itau = itau + 1
      if (itau > ntau_vcf) itau = itau - ntau_vcf

      forall(i = 1:nmmol)
         rOx(i,itau,:) = rcmO2(i,:)
         rN2(i,itau,:) = rcmN2(i,:)
         rCO2(i,itau,:) = rcmCO2(i,:)
      end forall

      if (icount <= ntau_vcf) RETURN
      isamp = isamp + 1
      do l = 1,ntau_vcf
         m = itau - l + 1
         if (m <= 0) m = m + ntau_vcf
         forall(i = 1:nmmol)
            rac1Ox(i,l) = rac1Ox(i,l) + rOx(i,itau,1)*rOx(i,m,1) &
                                      + rOx(i,itau,2)*rOx(i,m,2) &
                                      + rOx(i,itau,3)*rOx(i,m,3)
            rac1N2(i,l) = rac1N2(i,l) + rN2(i,itau,1)*rN2(i,m,1) &
                                      + rN2(i,itau,2)*rN2(i,m,2) &
                                      + rN2(i,itau,3)*rN2(i,m,3)
            rac1CO2(i,l) = rac1CO2(i,l) + rCO2(i,itau,1)*rCO2(i,m,1) &
                                        + rCO2(i,itau,2)*rCO2(i,m,2) &
                                        + rCO2(i,itau,3)*rCO2(i,m,3)
         end forall
      end do
   END SUBROUTINE vcf_samp

   SUBROUTINE vcf_print(nmmol,nfreq,dtfs)
      integer,intent(in):: nmmol,nfreq
      real(wp):: dtfs
      real(wp):: rnmal,tau
      integer:: j
      if (isamp <= 0) RETURN
      rnmal = 1.0/real(isamp,wp)
! sum individual molecule tracer VCFs
      rac1Ox_av(:) = 0.0
      rac1N2_av(:) = 0.0
      rac1CO2_av(:) = 0.0
      do j = 1,nmmol
         rac1Ox_av(:) = rac1Ox_av(:) + rac1Ox(j,:)
         rac1N2_av(:) = rac1N2_av(:) + rac1N2(j,:)
         rac1CO2_av(:) = rac1CO2_av(:) + rac1CO2(j,:)
      end do
!
      write(1004,"('# ',g15.7)") rac1Ox_av(1)*rnmal/nmmol
      write(1005,"('# ',g15.7)") rac1N2_av(1)*rnmal/nmmol
      write(1006,"('# ',g15.7)") rac1CO2_av(1)*rnmal/nmmol
      do j = 1, ntau_vcf
         tau = (j - 1)*nfreq*dtfs*0.001_wp
         write(1004,*) tau,rac1Ox_av(j)/rac1Ox_av(1)
         write(1005,*) tau,rac1N2_av(j)/rac1N2_av(1)
         write(1006,*) tau,rac1CO2_av(j)/rac1CO2_av(1)
      end do
      write(1004,'(/)')
      write(1005,'(/)')
      write(1006,'(/)')
   END SUBROUTINE vcf_print

END MODULE vcf_mod


PROGRAM SIMBOX
      USE commandline_mod
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
      USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE keating_parameters
      USE nlist_mod
      USE repul_energy_mod
      USE charges_mod
      USE force_en_el_mod
      USE insert_mod
      USE msd_mod
      USE vcf_mod
      implicit none
      integer,parameter:: natom_unit=24
      integer,parameter:: nunitc = 4
      integer,parameter:: number_moves = 1000000
      real(wp),parameter:: dtreal = 0.1_wp
      integer,parameter:: n_short = 4
      integer,parameter:: nequib = 100
      integer,parameter:: nsamp = 50
      integer,parameter:: nprint_all = 10000
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: mN = 14.00_wp
      real(wp),parameter:: mC = 12.00_wp
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,qj,sumqj,Temp_K
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp)::CL,time,Ek,EkO2,EkN2,EkCO2,rad_Ox(3),rad_N2(3),rad_CO2(3),rr(3),ULJEL,dt_short
      real(wp),allocatable:: massi(:,:),Ox_massi(:,:,:),N2_massi(:,:,:),CO2_massi(:,:,:)
      integer:: j,i,m,ii,ic = 0,natomc,istep,len,status,nat
      integer:: itmp,i_short,narg
      real(wp),allocatable:: coprxyz(:,:)
      real:: timep(0:10),sumtimep(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5
      character(len=132):: infile,carg
      logical:: restart
!
      ntau_msd = 10000
      ntau_vcf = 10000
!
      call cpu_time(timep(0))
      open(unit=24,file='final.pdb')
      open(unit=17,file='energy.out')
      open(unit=1001,file='msd_O2.out')
      open(unit=1002,file='msd_N2.out')
      open(unit=1003,file='msd_CO2.out')
      open(unit=1004,file='vcf_O2.out')
      open(unit=1005,file='vcf_N2.out')
      open(unit=1006,file='vcf_CO2.out')
!
      narg = command_argument_count()
      call get_command_argument (0, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', carg(1:len)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) carg(1:len),'  natom  atom_file  T[K]  {xyz_restart  vel_restart} '
         stop
      end if
!
      call get_command_argument (1, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 1
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) natom_max
      print *,'natom = ',natom_max
!
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 3
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) Temp_K
      call get_command_argument (4, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 4
         restart = .false.
      else
         restart = .true.
         write (*,*) 'restarting from ', infile(1:len)
         open(unit=34,file=trim(infile(1:len)))
         call get_command_argument (5, infile, len, status)
         if (status /= 0) then
            print *,'restart velocity file missing'
            stop
         else
            open(unit=35,file=trim(infile(1:len)))
         end if
      end if

!
!     CL = 0.7168_wp*0.162_wp/0.1552_wp
!     boxl = nunitc*CL
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!rcutel = boxl2
!rcutel2 = boxl2**2
      rad_Ox = (/sig_O_O2,sig_O_O2,0.0_wp/)*0.5_wp
      rad_N2 = (/sig_N_N2,sig_N_N2,0.0_wp/)*0.5_wp
      rad_CO2 = (/sig_C_CO2,sig_O_CO2,sig_O_CO2/)*0.5_wp
      n_O2 = 40
      n_N2 = 40
      n_CO2 = 40
      natom = natom_max
      T(1) = Temp_K
!
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
      nat = 2
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(Ox_xyz(1:n_O2,1:nat,1:3),N2_xyz(1:n_N2,1:nat,1:3),CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(r_O2_uc(n_O2,nat,3),r_N2_uc(n_n2,nat,3),r_CO2_uc(n_co2,3,3))
      allocate(dr_O2(n_O2,nat,3),dr_N2(n_n2,nat,3),dr_CO2(n_co2,3,3))
      allocate(rcm_O2_uc(n_O2,3),rcm_N2_uc(n_N2,3),rcm_CO2_uc(n_CO2,3))
      allocate(vcm_O2(n_O2,3),vcm_N2(n_N2,3),vcm_CO2(n_CO2,3))
      allocate(Imige_Ox_xyz(1:n_O2,1:3,1:3),Imige_N2_xyz(1:n_N2,1:3,1:3),Imige_CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(fxyz(1:natom_max,3),Ox_fxyz(1:n_O2,1:nat,1:3),N2_fxyz(1:n_N2,1:nat,1:3),CO2_fxyz(1:n_CO2,1:3,1:3))
      allocate(vxyz(1:natom_max,3),Ox_vxyz(1:n_O2,1:nat,1:3),N2_vxyz(1:n_N2,1:nat,1:3),CO2_vxyz(1:n_CO2,1:3,1:3))
      allocate(atom(1:natom_max),Ox_atom(1:n_O2,3),N2_atom(1:n_N2,3),CO2_atom(1:n_CO2,3))
      allocate(charge(natom),Ox_gas_charge(n_O2,3),N2_gas_charge(n_N2,3),CO2_gas_charge(n_CO2,3))
      allocate(proximity(natom_max,4))
      allocate(massi(natom_max,3),Ox_massi(1:n_O2,nat,3),N2_massi(1:n_N2,nat,3),CO2_massi(1:n_CO2,3,3))
      allocate(mol_vec(3,3))
      proximity = 0

      do i = 1,natom
         read (14,'(a6)',advance = 'no') c6
         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
!print *,i,atom(i),'   ',atom_name(atom(i))
         if (atom(i) == iSilicon) massi(i,:) = 1.0_wp/mSi
         if (atom(i) == iOxygen)  massi(i,:) = 1.0_wp/mOx
         if (atom(i) == iOxygenH)  massi(i,:) = 1.0_wp/mOH
      end do

      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
      do i = 1,natom
         read(14,'(a32)') ctmp
         do j = 1,ncmax(atom(i))
            c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
            read( unit=c5,fmt=* ) proximity(i,j)
         end do
      end do
      close(14)

      call assign_charge(1,natom)

write(*,*) 'sum charge = ',sum(charge(1:natom))
write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
write(*,*) 'n Si = ',count(atom == iSilicon)
write(*,*) 'n O  = ',count(atom == iOxygen)
write(*,*) 'n OH = ',count(atom == iOxygenH)
write(*,*) 'ntot = ',count(atom == iSilicon) + count(atom == iOxygen) + count(atom == iOxygenH)
!
! do i = 0,4
! write(888,*) i,qsi(i)
! end do
! do i = 1,natom
! write(888,*) i,atom_name(atom(i)),charge(i),noh(i)
! if (atom(i) ==iSilicon) then
!   sumqj = 0.0_wp
!   do j =1,4
!      ii = proximity(i,j)
!      if (ii == 0) STOP 'error in proximity'
!      qj = charge(ii)
!      if (atom(ii) ==iOxygen) qj = qj*0.5_wp
!      sumqj = sumqj + qj
!   end do
!   write(889,*) i,charge(i),sumqj,charge(i)+sumqj
! end if
! end do
!stop

      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom)
      call NEW_NLIST

!================== insertion procces==================================
      mol_vec(1,1) =  AOO*0.5_wp
      mol_vec(1,2) = -AOO*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_N2_O2(n_O2,Ox_xyz,1000000,rad_Ox)

      mol_vec(1,1) =  ANN*0.5_wp
      mol_vec(1,2) = -ANN*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_N2_O2(n_N2,N2_xyz,1000000,rad_N2)

      mol_vec(1,1) = 0.0_wp
      mol_vec(1,2) =  ACO2
      mol_vec(1,3) = -ACO2
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_CO2(n_CO2,CO2_xyz,1000000,rad_CO2)
!
      if (restart) then
         call read_box_frame(34,nmol = 40,print_all = .true.)
         call read_vel(35,nmol = 40,print_all = .true.)
      end if
!
! set up the uncorrected coords
      do i = 1, n_O2
         r_O2_uc(i,1,:) = Ox_xyz(i,1,:)
         rr = Ox_xyz(i,2,:) - Ox_xyz(i,1,:)
         call pbc(rr)
         r_O2_uc(i,2,:) = r_O2_uc(i,1,:) + rr
      end do
      do i = 1, n_n2
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
         call write_box_frame(imve,n_N2,print_all = .true.)
      end if
!========================================================================
!

      do i = 1,n_O2
        Ox_gas_charge(i,1) = q_O_O2
        Ox_gas_charge(i,2) = q_O_O2
        Ox_gas_charge(i,3) = -2.0_wp*q_O_O2
        Ox_atom(i,1) = iO_O2
        Ox_atom(i,2) = iO_O2
        Ox_massi(i,:,:) = 1.0_wp/mOx
      end do
      do i = 1,n_N2
        N2_gas_charge(i,1) = q_N_N2
        N2_gas_charge(i,2) = q_N_N2
        N2_gas_charge(i,3) = -2.0_wp*q_N_N2
        N2_atom(i,1) = iN_N2
        N2_atom(i,2) = iN_N2
        N2_massi(i,:,:) = 1.0_wp/mN
      end do
      do i = 1,n_CO2
        CO2_gas_charge(i,1) = q_C_CO2
        CO2_gas_charge(i,2) = q_O_CO2
        CO2_gas_charge(i,3) = q_O_CO2
        CO2_atom(i,1) = iC_CO2
        CO2_atom(i,2) = iO_CO2
        CO2_atom(i,3) = iO_CO2
        CO2_massi(i,1,:) = 1.0_wp/mC
        CO2_massi(i,2:3,:) = 1.0_wp/mOx
      end do

      call LJ_INIT

!     call comvel(natom,n_co2,n_o2,n_n2,t(1),vxyz,ox_vxyz,n2_vxyz,co2_vxyz)
!
!     vxyz = vxyz*sqrt(massi)
!     ox_vxyz = ox_vxyz*sqrt(ox_massi)
!     n2_vxyz = n2_vxyz*sqrt(n2_massi)
!     co2_vxyz = co2_vxyz*sqrt(co2_massi)
      vxyz = 0.0_wp
      Ox_vxyz = 0.0_wp
      N2_vxyz = 0.0_wp
      CO2_vxyz = 0.0_wp
      if (make_movie) then ! write out first frame of velocities
         iframe = 0
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_N2,print_all = .false.)
      end if

      call set_bond_angle_lists()
      call msd_init(40)
      call vcf_init(40)

      time = 0.0_wp
      dt = dtreal*dt_fs
      dt_short = dt/n_short     !0.25_wp*dt_fs

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      energy = 0.0_wp
!
!      call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!      call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!      call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! Electrostatic
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_n2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)
!      call FORCE_ENERGY_EL(1,natom,1,n_co2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)

!======================================================== main loop ===
      main_loop:do istep = 1,number_moves

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi
      Ox_vxyz = Ox_vxyz + 0.5_wp*dt*Ox_fxyz*Ox_massi
      N2_vxyz = N2_vxyz + 0.5_wp*dt*N2_fxyz*N2_massi
      CO2_vxyz = CO2_vxyz + 0.5_wp*dt*CO2_fxyz*CO2_massi

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      call cpu_time(timep(1))
      call force_keating
      call repul_energy2(1,natom)
      call O2_stretching
      call N2_stretching
      call CO2_stretching
      call CO2_angle_bending

!=============================== short loop ===
      short_loop: do i_short = 1, n_short
         coprxyz = rxyz
         call vverlet_a(dt_short,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)

         call gas_vverlet_a(dt_short,Ox_massi,Ox_xyz,Ox_vxyz,Ox_fxyz,dr_O2)
         call pbc(Ox_xyz)
         r_O2_uc = r_O2_uc + dr_O2

         call gas_vverlet_a(dt_short,N2_massi,N2_xyz,N2_vxyz,N2_fxyz,dr_N2)
         call pbc(N2_xyz)
         r_N2_uc = r_N2_uc + dr_N2

         call gas_vverlet_a(dt_short,CO2_massi,CO2_xyz,CO2_vxyz,CO2_fxyz,dr_CO2)
         call pbc(CO2_xyz)
         r_CO2_uc = r_CO2_uc + dr_CO2

         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3))) then
               !! print *,'NEW LIST'
               call NEW_NLIST
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         Ox_fxyz = 0.0_wp
         N2_fxyz = 0.0_wp
         CO2_fxyz = 0.0_wp
         energy = 0.0_wp
!call cpu_time(timep(4))
         call force_keating
!call cpu_time(timep(5)); sumtimep(4) = sumtimep(4) + timep(5)-timep(4)
         call repul_energy2(1,natom)
!call cpu_time(timep(6)); sumtimep(5) = sumtimep(5) + timep(6)-timep(5)
         call O2_stretching
         call N2_stretching
         call CO2_stretching
!call cpu_time(timep(7)); sumtimep(6) = sumtimep(6) + timep(7)-timep(6)
         call CO2_angle_bending
!call cpu_time(timep(8)); sumtimep(7) = sumtimep(7) + timep(8)-timep(7)
         call vverlet_b(dt_short,massi,vxyz,fxyz)
         call gas_vverlet_b(dt_short,Ox_massi,Ox_vxyz,Ox_fxyz)
         call gas_vverlet_b(dt_short,N2_massi,N2_vxyz,N2_fxyz)
         call gas_vverlet_b(dt_short,CO2_massi,CO2_vxyz,CO2_fxyz)
         !coprxyz = rxyz
!call cpu_time(timep(9)); sumtimep(8) = sumtimep(8) + timep(9)-timep(8)
!============================end short loop ===
      end do short_loop
      call cpu_time(timep(2)); sumtimep(1) = sumtimep(1) + timep(2) - timep(1)

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
!
!     call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!     call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!     call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! Electrostatic
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_n2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)
!      call FORCE_ENERGY_EL(1,natom,1,n_co2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)
!
      call cpu_time(timep(3)); sumtimep(2) = sumtimep(2) + timep(3) - timep(2)

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi
      Ox_vxyz = Ox_vxyz + 0.5_wp*dt*Ox_fxyz*Ox_massi
      N2_vxyz = N2_vxyz + 0.5_wp*dt*N2_fxyz*N2_massi
      CO2_vxyz = CO2_vxyz + 0.5_wp*dt*CO2_fxyz*CO2_massi

      time = time + dtreal
      if (istep <= nequib) then
         !  print *,sum(fxyz(1:natom,:),dim=1)
         Ek = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         do m = 1, n_O2
            do i = 1,2
               Ek = Ek + (0.5_wp/Ox_massi(m,i,1))*dot_product(Ox_vxyz(m,i,:),Ox_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_N2
            do i = 1,2
               Ek = Ek + (0.5_wp/N2_massi(m,i,1))*dot_product(N2_vxyz(m,i,:),N2_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_CO2
            do i = 1,3
               Ek = Ek + (0.5_wp/CO2_massi(m,i,1))*dot_product(CO2_vxyz(m,i,:),CO2_vxyz(m,i,:))
            end do
         end do

         T(2) = Ek*2.0_wp/( 3.0_wp*K_ev*(natom + 280) )
         if (mod(istep,nsamp) == 0) then
            write(*,*) istep,timep(2) - timep(1),timep(3) - timep(2),timep(3) - timep(0)
            write(17,'(8g18.8)') time, EK, energy + repulenergy*0.5_wp, Ek + energy + repulenergy*0.5_wp
            call flush(17)
         end if
         vxyz = vxyz*sqrt(T(1)/T(2))
         Ox_vxyz = Ox_vxyz*sqrt(T(1)/T(2))
         N2_vxyz = N2_vxyz*sqrt(T(1)/T(2))
         CO2_vxyz = CO2_vxyz*sqrt(T(1)/T(2))
      end if

      if ((mod(istep,nsamp) == 0) .AND. (istep >= nequib)) then
         forall(i = 1:n_O2)
            rcm_O2_uc(i,:) = ( r_O2_uc(i,1,:) + r_O2_uc(i,2,:) )*0.5_wp
            rcm_N2_uc(i,:) = ( r_N2_uc(i,1,:) + r_N2_uc(i,2,:) )*0.5_wp
            rcm_CO2_uc(i,:) = ( r_CO2_uc(i,1,:)*mC + r_CO2_uc(i,2,:)*mOx + r_CO2_uc(i,3,:)*mOx )/(mOx + mOx + mC)
            vcm_O2(i,:) = ( Ox_vxyz(i,1,:) + Ox_vxyz(i,2,:) )*0.5_wp
            vcm_N2(i,:) = ( N2_vxyz(i,1,:) + N2_vxyz(i,2,:) )*0.5_wp
            vcm_CO2(i,:) = ( CO2_vxyz(i,1,:)*mC + CO2_vxyz(i,2,:)*mOx + CO2_vxyz(i,3,:)*mOx )/(mOx + mOx + mC)
         end forall
         call msd_samp(n_O2,rcm_O2_uc,rcm_N2_uc,rcm_CO2_uc)
         call vcf_samp(n_O2,vcm_O2,vcm_N2,vcm_CO2)
         Ek = 0.0_wp
         EkO2 = 0.0_wp
         EkN2 = 0.0_wp
         EkCO2 = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         do m = 1, n_O2
            do i = 1,2
               EkO2 = EkO2 + (0.5_wp/Ox_massi(m,i,1))*dot_product(Ox_vxyz(m,i,:),Ox_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_N2
            do i = 1,2
               EkN2 = EkN2 + (0.5_wp/N2_massi(m,i,1))*dot_product(N2_vxyz(m,i,:),N2_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_CO2
            do i = 1,3
               EkCO2 = EkCO2 + (0.5_wp/CO2_massi(m,i,1))*dot_product(CO2_vxyz(m,i,:),CO2_vxyz(m,i,:))
            end do
         end do

         write (17,'(8g18.8)') time, Ek,EkO2,EkN2,EkCO2, Ek + EkO2 + EkN2 + EkCO2 + energy + repulenergy*0.5_wp
         call flush(17)
         write(*,*) istep,timep(2) - timep(1),timep(3) - timep(2),timep(3) - timep(0)
         if (make_movie) then
            iframe = iframe + 1
            call new_frame_file(imve,'frame',iframe)
            call write_box_frame(imve,n_n2,print_all = (mod(istep,nprint_all) == 0))
            call new_frame_file(imve,'frame_vel',iframe)
            call write_vel(imve,n_n2,print_all = (mod(istep,nprint_all) == 0))
         end if
      end if

!=====================================================end main loop ===
      end do main_loop
      call  msd_print(40,nsamp,dtreal)
      call  vcf_print(40,nsamp,dtreal)
      do i = 1,10
         print *,i,sumtimep(i)
      end do
      close(17,status='keep')
      close(1001,status='keep')
      close(1002,status='keep')
      close(1003,status='keep')
      close(1004,status='keep')
      close(1005,status='keep')
      close(1006,status='keep')
      call write_xyz(24,0)
      close(24,status='KEEP')
END PROGRAM simbox

