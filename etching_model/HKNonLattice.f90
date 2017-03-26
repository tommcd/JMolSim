
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
      integer:: i,ii,j,nnmax,iNodeNextL,inn,NodeLPmin,k
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

