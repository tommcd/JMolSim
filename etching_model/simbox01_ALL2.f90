
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
   integer,parameter,private :: ncmx = 2000
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



PROGRAM SIMBOX
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
      USE HKNonLattice_mod
      implicit none
      integer,parameter:: natomt = 1536
      real(wp):: r3(3)
      integer:: crossbond_x(natomt),crossbond_y(natomt),crossbond_z(natomt)
      real(wp):: CL,bl1,bl2
      integer:: j,i,m,n,ib,k,iat,ii,ic = 0
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      type(atom_list):: lst_positive_x,lst_negative_x,lstSi
      type(atom_list):: lst_positive_y,lst_negative_y
      type(atom_list):: lst_positive_z,lst_negative_z
      character(2):: ctmp
!
      open(unit=16,file='connect_struct.out')
      open(unit=15,file='neighbor.out')
      n_silicon_etch = 300
      natom = natomt
      natom_max = natomt
      nseed = 0
      allocate(rxyz(1:natom,3))
      allocate(atom(1:natom))
      allocate(proximity(natom,4))
      proximity = 0
      write (*,*) natom

! read input data
      open(unit=13,FILE = 'unitcell.dat')
      do i = 1,24
         read(13,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         if (i <= 8) then
            atom(i) = iSilicon
         else
            atom(i) = iOxygen
         end if
      end do
      close(13,status='keep')

!------------------------------------------------------------
      CL = 40.0_wp
      boxl = 4.0_wp*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      bl1 = CL*sqrt(2.0_wp)/8.0_wp
      bl2 = CL*sqrt(1.5_wp)/4.0_wp

!
! Creating the simulation box
!
      open(unit=14,file='simbox.out')
      do k = 1,4
      do j = 1,4
      do i = 1,4
      do m = 1,24
         write(14,*) atom_name(atom(m)),rxyz(m,1) + (i - 1)*CL, &
                     rxyz(m,2) + (j - 1)*CL,rxyz(m,3) + (k - 1)*CL
      end do
      end do
      end do
      end do
      rewind(14)
! origin at centre of box
      do i = 1,natom
         READ(14,*) ctmp, rxyz(i,1), rxyz(i,2), rxyz(i,3)
         atom(i) = name2atom(ctmp)
         rxyz(i,1) = rxyz(i,1) - 1.5_wp*CL
         rxyz(i,2) = rxyz(i,2) - 1.5_wp*CL
         rxyz(i,3) = rxyz(i,3) - 1.5_wp*CL
!--------Periodic boundary conditions
         rxyz(i,1) = rxyz(i,1) - boxl*anint(rxyz(i,1)*boxli)
         rxyz(i,2) = rxyz(i,2) - boxl*anint(rxyz(i,2)*boxli)
         rxyz(i,3) = rxyz(i,3) - boxl*anint(rxyz(i,3)*boxli)
      end do
      close(14,status='KEEP')
call write_xyz(991,0)

! determine connectivity
      do k = 1,natom
         n = 0
         do j = 1,natom
            if (j == k) cycle
            if (atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
                n = n + 1
                proximity(k,n) = j
            end if
         end do
      end do
      do i = 1,natom
         write (15,*) atom_name(atom(i)),proximity(i,1:4)
      end do
! count the clusters
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
         stop 'there should be one cluster'
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

! Make array of silicon atoms
      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == iSilicon) call add_to_list(lstSi,i)
      end do

!
! Perfom the etching procedure
!
      do i = 1, n_silicon_etch

         ! pick a Silicon to remove
         !call rand_from_list(lstSi,iat)
         do
            iat = int(rand()*natom) + 1
            if (atom(iat) == iSilicon) exit
         end do

         call delete_atom(iat)

         !call remove_from_list(lstSi,iat)
         !if (atom(iat) == iSilicon) then
         !   where(lstSi%i == (natom + 1)) lstSi%i = iat
         !end if
         ! check for errors
         if (ANY(proximity(1:natom,:) > natom)) then
            print *,'error: proximity(1:natom,:) > natom ',i
            stop
         end if
         if (ANY(proximity(natom + 1:,:) /= 0)) then
            print *,'error: proximity(natom + 1:,:) /= 0 ',i
            stop
         end if
      end do

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if

!
! Make list of bonds
      call bond_list

!
! Make arrays of crossing bonds
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
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

!
! Make arrays of interface atoms
!
      write(*,*) n_crossbond_x
      write(*,*) n_crossbond_y
      write(*,*) n_crossbond_z
      if ((n_crossbond_x == 0).or. &
         (n_crossbond_y == 0).or. &
         (n_crossbond_z == 0)) then
         print *,'try again'
         stop
      end if

!
      do i = 1,n_crossbond_x
         ib = crossbond_x(i)
         if (rxyz(ibond(1,ib),1) > 0.0) then
            call add_to_list(lst_positive_x,ibond(1,ib))
            call add_to_list(lst_negative_x,ibond(2,ib))
            if (rxyz(ibond(2,ib),1) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),1),rxyz(ibond(1,ib),1)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_x,ibond(1,ib))
            call add_to_list(lst_positive_x,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_y
         ib = crossbond_y(i)
         if (rxyz(ibond(1,ib),2) > 0.0) then
            call add_to_list(lst_positive_y,ibond(1,ib))
            call add_to_list(lst_negative_y,ibond(2,ib))
            if (rxyz(ibond(2,ib),2) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),2),rxyz(ibond(1,ib),2)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_y,ibond(1,ib))
            call add_to_list(lst_positive_y,ibond(2,ib))
         end if
      end do
!
      do i = 1,n_crossbond_z
         ib = crossbond_z(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_positive_z,ibond(1,ib))
            call add_to_list(lst_negative_z,ibond(2,ib))
            if (rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_z,ibond(1,ib))
            call add_to_list(lst_positive_z,ibond(2,ib))
         end if
      end do
!

!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

!
! Finding a cluster which crosses all faces
!
      ic = 0
      do i = 1,n_crossbond_x

         connected_x = .false.
         do j = 1,n_crossbond_x
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_x%i(j))) then
               connected_x = .true.
               exit
            end if
         end do
         if (.not.connected_x) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.not.connected_y) CYCLE

         connected_y = .false.
         do j = 1,n_crossbond_y
            if (atomL(lst_positive_x%i(i)) == atomL(lst_positive_y%i(j))) then
               connected_y = .true.
               exit
            end if
         end do
         if (.not.connected_y) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(i)) == atomL(lst_negative_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.not.connected_z) CYCLE

         connected_z = .false.
         do j = 1,n_crossbond_z
            if (atomL(lst_positive_x%i(i)) == atomL(lst_positive_z%i(j))) then
               connected_z = .true.
               cycle
            end if
         end do
         if (.not.connected_z) CYCLE

         spanning_cluster = connected_x.AND.connected_y.AND.connected_z
         if (spanning_cluster) then
            ic = atomL(lst_positive_x%i(i))
            EXIT
         end if
      end do
      if (.not.spanning_cluster) then
         stop 'no spanning cluster found'
      end if

      if (spanning_cluster) then
      do j = 1,natom
         if (atomL(j) == ic) then
            write(16,*) atom_name(atom(j)),rxyz(j,1),rxyz(j,2),rxyz(j,3)
         end if
      end do
      end if
      CLOSE(UNIT = 16,STATUS = 'KEEP')

CONTAINS

   pure FUNCTION len_3d(r)
      real(wp):: len_3d
      real(wp),intent(in):: r(3)
        len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
   END FUNCTION

END PROGRAM SIMBOX

