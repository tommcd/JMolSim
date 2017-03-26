
include 'HKNonLattice2.f90'
include 'qsort.f90'


program test_HKNonLattice
      USE HKNonLattice_mod
      implicit none
      integer,allocatable:: Node(:),NodeNext(:,:),NodeL(:),atom(:)
!!    integer,allocatable:: LinkS(:),LinksOfNode(:,:),NodeS(:)
      integer:: fid,N_Node,N_Link,i,j,NumberOfClusters,it
      character(132):: c2*2,ctmp*32,fin0
      integer,parameter:: ncmax(0:1) = (/2,4/)
      integer,external:: iargc
!
      fid = 10
!     if (command_argument_count() /= 1) then
      if (iargc() /= 1) then
        write(*,*)'1 argument(s) required :'
        write(*,*)'  1. input file'
        stop
      end if

!     call get_command_argument(1,fin0)   ;write(*,*) trim(fin0)
      call getarg(1,fin0)   ;write(*,*) trim(fin0)
      open(unit=fid,file=trim(fin0))
      read(fid,*) N_Node   ;print *,'N_Node = ',N_Node
      allocate( Node(N_Node) )
      allocate( NodeNext(N_Node,4) )
      allocate( NodeL(N_Node))
      allocate( atom(N_Node))
!
!!    allocate( LinksOfNode(N_Node,4) )
!!    allocate( NodeS(N_Node) )
!
      NodeNext = 0
      do i = 1,N_Node
         read(fid,*) atom(i),it
!!       read(fid,*) Node(i), NodeS(i), NodeNext(i,1:4)
!!       read(fid,*) Node(i), NodeNext(i,1:4)
      end do
      do i = 1,N_Node
         read(fid,*) it,ctmp,it,(NodeNext(i,j),j = 1,ncmax(atom(i)))
      end do

!!      read(fid,*) N_Link   ;print *,'N_Link = ',N_Link
!!      allocate( LinkS(N_Link) )
!!      do i = 1, N_Link
!!         read(fid,*) ii,LinkS(i)
!!      end do
!!      read(fid,*) N_Node   ;print *,'N_Node = ',N_Node
!!      do i = 1, N_Node
!!         read(fid,*) ii, LinksOfNode(ii,1:4)
!!      end do
!
print *,'N_Node = ',N_Node
      write(77,*) 'NodeNext = '
      do i = 1,N_Node
         write(77,'(6i6)') i,NodeNext(i,:)
      end do
stop
      call HKNonLattice(N_Node,NodeNext,NumberOfClusters,NodeL)
      write(77,*) 'NumberOfClusters = ',NumberOfClusters
      write(77,*) 'NodeL: ',NodeL
!
end program test_HKNonLattice

