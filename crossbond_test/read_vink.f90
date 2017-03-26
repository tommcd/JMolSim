
include 'precision_mod.f90'
include 'global_vars.f90'
include 'sort_mod.f90'
include 'files_mod.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'nlist_mod.f90'
include 'bond_list_mod.f90'
include 'frames_mod.f90'
include 'readline_mod.f90'
include 'rand_mod.f90'
include 'matvec3d_mod.f90'
include 'HKNonLattice2.f90'


PROGRAM crossbond_correct
      USE precision_mod, only: wp
      USE global_vars
      USE bond_list_mod
      USE constants_mod
      USE atom_types
      USE connectivity_mod
      USE coordinates_mod
      USE files_mod
      USE seaton_mod
      USE readline_mod
      USE nlist_mod
      USE frames_mod
      USE HKNonLattice_mod
      IMPLICIT NONE
      integer,parameter:: nl = 1 
      integer,parameter:: nunitc(3) = 2*nl + 1
      integer:: ip,itmp,ifirst,nsi,nox
      integer:: i,ii,j,jj,k,kk,m,n
      integer:: iu,ios,ia,ig,ilast,iatom,natom_orig,ix,iy,iz,ib
      integer:: n_crossbond,ncross,nc,ic,neigh2,neighmx=125
      integer,allocatable:: atomc(:),proximity_temp(:,:)
      integer,allocatable:: cross_atom(:),atom_temp(:),crossbond(:)
      real(wp):: timep(0:10),bond_lens(4)
      real(wp):: r3(3),rj(3),dr(3),bondl_orig=1.55_wp
      real(wp),allocatable:: rxyz_temp(:,:)
      real(wp):: CL,CL2,x,y,z,boxl_0(3),en
      character(len=32):: carg,word,ctmp
      character(6):: c6
      character(5):: c5
      logical:: consistent
!
      CALL get_free_file_unit(iu)
      open(unit=iu,file='60k.sio2')
      read(iu,*) CL2  ; print *,'BOXL2 = ',CL2
      CL = 2.0_wp*CL2  ; print *,'BOXL  = ',CL
      read(iu,*) en  ; print *,'en = ',en
      read(iu,*) nsi
      read(iu,*) nox
      natom = nsi + nox
      natom_max = natom
      allocate(rxyz(natom_max,3),proximity(natom_max,4))
      allocate(atom(natom_max))
      do i = 1,nsi
         read(iu,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         atom(i) = iSilicon
      end do
      do i = nsi+1,nsi+nox
         read(iu,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         atom(i) = iOxygen
      end do
!
      boxl = CL
      boxli = 1.0_wp/boxl
      boxl_0(:) = boxl(:)
      boxl2 = boxl/2.0_wp
!

      CALL cpu_time(timep(0))
!
!     CALL get_free_file_unit(iu)
!     do i=1,natom
!        read (iu,*) c6,itmp,ctmp,itmp,rxyz_temp(i,:)
!        atom_temp(i) = name2atom(trim(ctmp))
!     end do
!     do i=1,natom
!        read(iu,'(a32)')ctmp
!        do j = 1,ncmax(atom_temp(i))
!           c5 = ctmp(6+5*(j)+1:6+5*(j)+5)
!           !print*, c5
!           read( unit=c5,fmt=* ) proximity_temp(i,j)
!        end do
!     end do
!     CLOSE(iu, STATUS='KEEP')
      
      proximity = 0   
      call INIT_NLIST(boxl,5.0_wp)
      call NEW_NLIST
!
! Loop over the Oxygen atoms
      oxygen_loop: do i = nsi+1, nsi+nox
         r3(1:3) = rxyz(i,1:3)
         ic = CELL(r3)
         CALL NEIGCELL(ic,2,neigh,ncell)
         n = 0
         cell_loop: do jj = 1, neigh
            nc = ncell(jj)
            if(nc==0) cycle cell_loop
            j=HOC(nc)
            cell_atom_loop: do while (j/=0)   ! loop over cell Silicon atoms
               if(atom(j) /= iSilicon) GOTO 100
               if (j==i) GOTO 100
               dr = r3-rxyz(j,1:3)
               CALL pbc(dr)
               if (dot_product(dr,dr) < (1.1*bondl_SiO)**2) then 
                  n=n+1
                  proximity(i,n) = j
                  call set_proximity(j,0,i)
                  if (n > 2) then   
                     print*, 'defect zone,n = ', n
                     print*, 'prox(',i,',:) = ', proximity(i,:)
                     stop
                  end if
               end if
100            CONTINUE
               j=LL(j)
            end do cell_atom_loop
         end do cell_loop 
      end do oxygen_loop
      
      CALL cpu_time(timep(1))
      print*, 'time taken = ', timep(1) - timep(0)
      
      call check_proximity(consistent,ifirst,k)
      if (.NOT. consistent) then
         print*, 'not consistent'
         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
         print*, 'proximity(',k,',:)=', proximity(k,:)
      end if

      do i = 1,nsi
         if(count(proximity(i,:)>0) /= 4)then
            print *,'proximity error'
            print *,'i = ',i
            print *,(proximity(i,j),j=1,4)
         end if
      end do
      do i = nsi+1,nsi+nox
         if(count(proximity(i,:)>0) /= 2)then
            print *,'proximity error'
            print *,'i = ',i
            print *,(proximity(i,j),j=1,4)
         end if
         if(is_2Si_ring(i))then
            print *,i,' is in a 2-Si ring'
         end if
      end do

      call check_proximity2(consistent,ifirst)
      if (.NOT. consistent) then
         print*, 'not consistent 2'
         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
      end if


      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if(n_cluster /= 1)then
         stop 'there should be one cluster'
      end if
      CALL get_free_file_unit(iu)
      CALL new_frame_file(iu,'60k_SiO2',0)
      CALL write_frame(iu,1,natom)
      close(iu)
END PROGRAM crossbond_correct

