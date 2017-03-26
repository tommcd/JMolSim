
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



PROGRAM perc2mol
      USE HKNonLattice_mod
      implicit none
      integer,parameter:: wp = kind(1.0d0)
      integer,parameter:: N = 14
      integer,parameter:: nbin = N**3
      real(wp),parameter:: angstrom = 10.0_wp
      integer,parameter:: full = 0,free = 1
      real(wp):: Kcrit = 0.035_wp
      integer:: i,j,k,narg,len,status,nb,iu
      integer:: ip,im,jp,jm,kp,km,na
      real(wp):: x,y,z
      real(wp):: henryk(N,N,N)
      real(wp):: rxyz(nbin,3),hk(nbin)
      integer:: a(N,N,N) = full
      integer:: b(N,N,N) = 0
      integer:: nc(nbin)=0,proximity(nbin,6)=0
      character(len=132) infile,outfile,ctmp
      real(wp):: boxl,boxl2 !,boxli
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
!     boxli = 1.0_wp/boxl
      call Init_HKNonLattice(nbin)
!
      narg = command_argument_count()
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      if(narg < 3)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' Kcrit infile outfile'
         stop
      end if
!
      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) Kcrit ; print *,'Kcrit = ',Kcrit
!
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, outfile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', outfile(1:len)
      iu = 15
      open(unit=iu,file=trim(outfile(1:len)))

!
      na = 0
      do k = 1,N
         do j = 1,N
            do i = 1,N
               read (14,*) x,y,z,henryk(i,j,k)
!print *,x,y,z,henryk(i,j,k)
               if (henryk(i,j,k) > Kcrit) then
                  na = na + 1
                  a(i,j,k) = free
                  rxyz(na,:) = (/ x,y,z /) - (/ boxl2,boxl2,boxl2 /)
                  hk(na) = henryk(i,j,k)/Kcrit
               end if
            end do
         end do
      end do
!
      do k = 1,N
         do j = 1,N               
            write(*,'(20i1)')(a(i,j,k),i=1,N)
         end do
         write(*,*)
      end do
      write(*,*) 'bins with K > Kcrit ',count(a == free)

      b = 0
      na = 0
      do k = 1,N
      do j = 1,N
      do i = 1,N
         ! ib = i + (j-1)*N + (k-1)*N*N
         if (a(i,j,k) == full) CYCLE
         na = na + 1
         b(i,j,k) = na
      end do         
      end do
      end do
!
      na = 0
      do k = 1,N
         kp = k + 1
         if (kp > N) kp = 1
         km = k - 1
         if (km < 1) km = N
         do j = 1,N
            jp = j + 1
            if (jp > N) jp = 1
            jm = j - 1
            if (jm < 1) jm = N
            do i = 1,N
               if (a(i,j,k) == full) CYCLE
               na = na + 1
               ip = i + 1
               if (ip > N) ip = 1
               im = i - 1
               if (im < 1) im = N
               ! connectivity
               nb = 0
               if (b(ip,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(ip,j,k)
               end if
               if (b(im,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(im,j,k)
               end if
               if (b(i,jp,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jp,k)
               end if
               if (b(i,jm,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jm,k)
               end if
               if (b(i,j,kp) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,j,kp)
               end if
               if (b(i,j,km) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,j,km)
               end if
               nc(na) = nb
            end do         
         end do
      end do
!
      call HKNonLattice(na,proximity,n_cluster,atomL)
      call analyse_cluster(na)
      call print_cluster(na,6)
!
      do i = 1,na
         write(iu,110)'HETATM',i,' C  ','   ',i,rxyz(i,1:3)*angstrom,hk(i),hk(i)
      end do
      do i = 1,na
         if ( all(proximity(i,:) == 0) ) cycle
         write(iu,"(a6,7i5)")'CONECT',i,(proximity(i,j),j = 1,nc(i))
      end do
      write(iu,110)'HETATM',99991,'AR','   ',99991,-boxl2*angstrom,-boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99992,'AR','   ',99992,-boxl2*angstrom, boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2*angstrom,-boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2*angstrom, boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99995,'AR','   ',99995,-boxl2*angstrom,-boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99996,'AR','   ',99996,-boxl2*angstrom, boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2*angstrom,-boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2*angstrom, boxl2*angstrom, boxl2*angstrom
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
110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3,2f6.2)
contains

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
end program

