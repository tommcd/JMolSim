INCLUDE 'precision_mod.f90'
INCLUDE 'math_const_mod.f90'
INCLUDE 'phys_const_mod.f90'
INCLUDE 'global_vars_mod_KMC.f90'
INCLUDE 'files_mod.f90'
INCLUDE 'sort_mod.f90'
INCLUDE 'atom_types_mod.f90'
INCLUDE 'coordinates_mod.f90'
INCLUDE 'connectivity_mod.f90'
INCLUDE 'nlist_mod3.f90'
INCLUDE 'HKNonLattice2.f90'
INCLUDE 'check_structure_mod_Silica.f90'

!>>INCLUDE 'readline_mod.f90'
MODULE readline_mod
    IMPLICIT NONE
    save
    integer,parameter :: maxlnt = 132
    character(len=*),parameter :: fmtt = '(a132)'
    character(len=maxlnt) :: line
    integer :: istart,lnt
    public:: getword,getline,line,istart,lnt
!
CONTAINS

  SUBROUTINE getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
  END SUBROUTINE getline

  SUBROUTINE getword(word,ios)
      character(len=*),intent(out) :: word
      integer,intent(out):: ios
      integer :: i,iend,i1
      ios = 0
      word(:)=''
      if(istart > lnt)then
         ios = -1
         return
      end if
      i1 = 0
      do i = istart,lnt
         if(line(i:i) /= ' ')then
            i1 = i
            exit
         end if
      end do
      if(i1 == 0) then
         ios = -2
         return
      end if
      iend = 1
      do i = i1+1,lnt
         if(line(i:i) == ' ' .or. i == lnt)then
            iend = i
            exit
         end if
      end do
      if (i1 >= lnt) iend = i1
      word(1:iend-i1+1) = line(i1:iend)
      istart = i+1
      return
  END SUBROUTINE getword
END MODULE readline_mod

INCLUDE 'matvec3d_mod.f90'

!>> INCLUDE 'pdb_read_mod.f90'
MODULE pdb_read_mod
   USE precision_mod
   USE global_vars_mod
   USE files_mod
   USE atom_types_mod
   USE coordinates_mod
   USE connectivity_mod
   IMPLICIT NONE
   SAVE
   integer,parameter :: maxlnt = 132
   character(len=*),parameter :: fmtt = '(a132)'
   character(len=maxlnt) :: line
   integer :: istart,lnt
   public:: getword,getline,line,istart,lnt
CONTAINS

   SUBROUTINE pdb_read(iu,natom)
      integer:: iargc,iatom,j
      integer:: ia,ig,iu,io,i,ios
      integer,intent(out):: natom
      character(len=132):: fin0,fout0,fout1
      character(len=80):: word,carg
      character(len=32):: c5,atomname
      natom = 0
      iatom = 0
      call getline(iu,ios)
      call getword(carg,ios)
      call getword(carg,ios)
      read(unit=carg,fmt=*) natom
      call getline(iu,ios)
      call getword(carg,ios)
      call getword(carg,ios)
      read(unit=carg,fmt=*) boxl(1)
      boxl(1:3) = boxl(1)
      allocate(rxyz(natom,3),proximity(natom,4),atom(natom))
      rxyz = 0.0_wp
      proximity = 0
      i=0
      lines: do
         i=i+1
         call getline(iu,ios)
         if(ios /= 0) exit
         call getword(word,ios)
         word = line(1:6)
         select case(trim(word))
         case('HETATM')
            if (i<10000) call getword(carg,ios)
            carg = line(6+1:6+5)
            read(unit=carg,fmt=*)ia
            call getword(carg,ios)
!            read(unit=carg,fmt=*)atom(name2atom(atom_name(i))
            read(unit=carg,fmt=*)atomname
            atom(i) = name2atom(trim(atomname))
            call getword(carg,ios)
            read(unit=carg,fmt=*)ig
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,1)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,2)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(i,3)
        case('CONECT')
            call getword(carg,ios)
            c5 = line(6+1:6+5)
            read(unit=c5,fmt=*)iatom
            do j=1,ncmax(atom(iatom))
               c5 = line(6+5*(j)+1:6+5*(j)+5)
               read(unit=c5,fmt=*) proximity(iatom,j)
            end do
            if(iatom==natom) EXIT lines
        case('END')
            EXIT
        case default
            write(*,*) 'error : ',trim(word),' is not on list'
        end select
      end do lines
   CLOSE(iu)

   END SUBROUTINE pdb_read

   SUBROUTINE getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
   END SUBROUTINE getline

   SUBROUTINE getword(word,ios)
      character(len=*),intent(out) :: word
      integer,intent(out):: ios
      integer :: i,iend,i1
      ios = 0
      word(:)=''
      if(istart > lnt)then
         ios = -1
         return
      end if
      i1 = 0
      do i = istart,lnt
         if(line(i:i) /= '')then
            i1 = i
            exit
         end if
      end do
      if(i1 == 0) then
         ios = -2
         return
      end if
      iend = 1
      do i = i1+1,lnt
         if(line(i:i) == '' .or. i == lnt)then
            iend = i
            exit
         end if
      end do
      if (i1 >= lnt) iend = i1
      word(1:iend-i1+1) = line(i1:iend)
      istart = i+1
      return
   END SUBROUTINE  getword

END MODULE pdb_read_mod

!>>INCLUDE 'bond_list_mod.f90'

MODULE bond_list_mod
   use precision_mod
   implicit none
   integer,allocatable:: ibond(:,:),nbond(:)
   integer,allocatable:: iang(:,:)
   real(wp),allocatable:: kang(:),ctheta(:)
   integer:: nbondtot,nang
CONTAINS
!
   SUBROUTINE bond_list
      USE global_vars_mod, only: natom,natom_max
      USE connectivity_mod, only: proximity
      USE atom_types_mod
!      USE Keating_parameters
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

!>>INCLUDE 'crossbond_mod.f90'
MODULE crossbond_mod

!    READ IN THE INITIAL DATA
!    WHEN HAVE THE DATA THEN WE NEED TO IDENTIFY ALL ATOMS INVOLVED IN CROSSBONDS
!    THIS WILL BE DONE WITH THE CROSSBOND CODE IN DEVELOPEMENT. IT IS NOW NOT
!    REQUIRED TO ACTUALLY SEGREGATE THE BONDS INTO THE DIFFERENT DIRECTIONS AS THE
!    PROXIMITY VALUES OF EACH WILL NOW BE SET TO ZERO AND CALCULATED USING THE
!    NEIGHBOUR LIST ALGORITHM.
!    ALTERNATIVELY TO THE ORIGINAL APPROACH TO THE PROBLEM OF CORRECTING THE
!    CROSSBONDS TO THEIR IMAGES IN THE VARIOUS DIRECTIONS, THIS WILL IMAGE ALL ATOM
!    COORDINATES AND ALL RESPECTIVE PROXIMITY VALUES EXCEPT THOSE OF ATOMS INVOLVED
!    IN CROSSBONDS. IN THE LATTER CASE THE PROXIMITY VALUES OF THE ATOM ARE SET TO
!    ZERO.
!    FOLLOWING THIS THE SITUATION EXISTS THAT ALL THE CONNECTING ATOMS HAVE NO
!    CONNECTIVITY. THESE ARE LOOPED OVER AND AN ATTEMPT IS MADE TO CORRECT THE
!    CONNECTIVITY OF THE SYSTEM. THIS IS DONE USING THE INITIAL METHOD OF ATOMS
!    BEING WITHIN AN SI-O BOND LENGTH AWAY FROM EACH OTHER. IF IT HAPPENS THAT ONE
!    OF THESE ATOMS HAS A NUMBER OF NEIGHBOURS CLOSER THAN THE CUT-OFF DISTANCE
!    WHICH IS HIGHER THAN THE DESIRED COORDINATION THEN WE GO BACK TO THE SAME ATOM
!    IN ORIGINAL BOX AND TAKE THE "EXACT" DISTANCES BETWEEN IT AND IT'S NEIGHBOURS.
!    USING THESE DISTANCES WE SEARCH FOR NEIGHBOURS WITH CORRESPONDING EXACT
!    DISTANCES FROM TARGET ATOM.


!    if imaging a box, then it is required that the crossbond mod should be
!    initialised before the imaging takes place. then the crossbond correction
!    subroutine should be called after the imaging of the coordinates and
!    connectivity takes place.
!    NOTE:- the initial box must be a cube for this to work

   USE precision_mod
   USE coordinates_mod
   USE connectivity_mod
   IMPLICIT NONE
   integer:: ncross
   integer,allocatable:: proximity_temp(:,:)
   integer,allocatable:: cross_atom(:)
CONTAINS

   SUBROUTINE init_crossbond(n_image,boxl)
      USE bond_list_mod
      USE nlist_mod
      integer,intent(in):: n_image !total number of repeated original cells
      real(wp),intent(in):: boxl(3)

      integer:: i,j,n_crossbond
      integer,allocatable:: crossbond(:),cross_atom(:)
      real(wp):: r3(3)

      ! this will set the various values of initial box length and the proximity
      ! values required for the initial deletion of the values
      ! detail on the input box's size
      boxl2 = boxl/2.0_wp
      ! CHECK FOR CROSSBONDS.
      call bond_list
      if (.not.allocated(crossbond)) allocate(crossbond(nbondtot))
      crossbond = 0
      n_crossbond = 0
      do i=1,nbondtot ! nbondtot calc'ed in bond_list
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         r3 = SQRT(dot_product(r3,r3))
         if (  r3(1) > boxl2(1) .OR. r3(2) > boxl2(2) .OR. &
               r3(3) > boxl2(3)  ) then
            n_crossbond = n_crossbond + 1
            crossbond(n_crossbond) = i
         end if
      end do
!    store all of the atom numbers that are involved in crossbonds and the
!    connectivity of these atoms
      allocate(cross_atom((n_image)*nbondtot*2))
      j=1
      ncross = 0
      do i=1,n_crossbond
         cross_atom(j) = ibond(1,crossbond(i))
         cross_atom(j+1) = ibond(2,crossbond(i))
         ncross = ncross + 2
         j=j+2
      end do

      proximity_temp(1:natom,:) = proximity(1:natom,:)

      do i =1, ncross
         proximity(cross_atom(i),:) = 0
      end do
   END SUBROUTINE init_crossbond


   SUBROUTINE crossbond_correct(nix,niy,niz,boxl_0,natom_orig)
      !!! NOTE:- only works for atoms with 4 or less neighbours 
      USE precision_mod, only: wp
      USE global_vars_mod
      USE bond_list_mod
      USE atom_types_mod
      USE connectivity_mod
      USE files_mod
      USE nlist_mod
      ! nix, niy, & niz are the number of images in each respective direction so
      ! need to have this passed in...probably the easiest thing to do
      integer,intent(in):: natom_orig
      real(wp),intent(in):: nix,niy,niz,boxl_0(3)

      integer:: ic,nc,iatom
      integer:: i,ii,k,kk,m,n !counters and dummy args
      integer:: ifirst,neigh,ncell(125)
      real(wp):: r3(3),bond_lens(4),dr(3)
      logical:: consistent

      !  need to specify the changes in the box length directions for
      !  initialisation of the neighbour list.
      boxl(1:3) = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
      boxl2(1:3) = boxl(1:3)/2.0_wp

      call  DELETE_NLIST()     

      print*, 'nix = ',nix
      print*, 'niy = ',niy
      print*, 'niz = ',niz

      call INIT_NLIST(-boxl2,boxl2,5.0_wp,natom)
      call NEW_NLIST(1,natom)

      ! the following is to get the connectivity of the edge atoms using the
      ! adjusted nlist
      main_loop: do i = 1, natom
         if ( .NOT. ALL( proximity(i,:)==0 ) ) cycle
         r3(1:3) = rxyz(i,1:3)
         ic = CELL(r3)
         CALL NEIGCELL(ic,2,neigh,ncell)
         n=0
         if (i>natom_orig) then
            iatom = mod(i,natom_orig)
         elseif (i<=natom_orig) then
            iatom = i
         end if
         ! need to have an array containing the connectivity
         ! of the target atom i. iatom is the original of
         ! the image under investigation.
         boxl = (/boxl_0(1),boxl_0(2),boxl_0(3)/)
         boxl2 = boxl/2.0_wp
         bond_lens = 0.0_wp
         do k = 1, ncmax(atom(iatom))
            dr = rxyz(iatom,:) - rxyz(proximity_temp(iatom,k),:)
            CALL pbc(dr)
            bond_lens(k) = dot_product(dr,dr)
         end do
         !reset boxl to proper value
         boxl = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
         boxl2 = boxl/2.0_wp
            inner_cell_loop: do kk=1,neigh
               nc=ncell(kk)
               if (nc==0) cycle inner_cell_loop
               m=HOC(nc)
               inner_cell_atom_loop: do while (m/=0)
                  dr = r3-rxyz(m,1:3)
                  if (m==i) GOTO 20
                  CALL pbc(dr)
                  do ii = 1, ncmax(atom(i))
                     if (ABS( dot_product(dr,dr)-bond_lens(ii) )<=0.00001_wp) EXIT
                     if (ii==ncmax(atom(i))) GOTO 20
                  end do
                  n=n+1
                  proximity(i,n) = m
                  if (n>4) stop 'need better condition'
                  if (n==ncmax(atom(i))) cycle main_loop
20                CONTINUE
                  m=LL(m)
               end do inner_cell_atom_loop
            end do inner_cell_loop
      end do main_loop

      boxl = (/boxl_0(1)*nix,boxl_0(2)*niy,boxl_0(3)*niz/)
      boxl2 = boxl/2.0_wp
      call bond_list

! turn on the following checks if you don't trust the correction
!      call check_proximity(consistent,ifirst,k)
!      if (.NOT. consistent) then
!         print*, ''
!         print*, 'not consistent'
!         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
!         print*, 'proximity(',k,',:)=', proximity(k,:)
!      end if
!
!      call check_proximity2(consistent,ifirst,k)
!      if (.NOT. consistent) then
!         print*, ''
!         print*, 'not consistent_2'
!         print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
!         print*, 'atom_name(atom(ifirst)) = ', atom_name(atom(ifirst))
!         print*, 'atom_name(atom(k)) = ', atom_name(atom(k))
!      end if
   END SUBROUTINE crossbond_correct
END MODULE crossbond_mod


INCLUDE 'command_line_mod.f90'


PROGRAM etch_amorphous
      USE precision_mod
      USE atom_types_mod
      USE files_mod
      USE coordinates_mod
      USE connectivity_mod
      USE global_vars_mod
      USE nlist_mod
      USE HKNonlattice_mod
      USE readline_mod
      USE sort_mod
      USE matvec3d_mod
      USE pdb_read_mod
      USE crossbond_mod
      USE COMMAND_LINE_MOD
      USE check_structure_mod

      IMPLICIT NONE
      integer,parameter:: iBorate=0
      integer,parameter:: iSilicate=1
      real(wp),parameter:: bondl_SiO = 1.60_wp
      character(2),parameter:: atom_name2(0:1)=(/ 'B ','Si' /)
   
      integer:: narg,leng,stat
      integer:: i,j,k,ii,jj,iu,ic,nc
      integer:: floater
      integer:: ia,ig,ios,iatom
      integer:: ps_natom,ndel,tmp
      integer:: n_OH,coord_type,nat_temp
      integer:: ifirst,i2,ucoord1,ucoord2,ucoord3,MaxCluster
      integer:: nSi_1_OH,nSi_1_OH2,nSi_1_OH3
      integer:: nSi_2_OH,nSi_2_OH2,nSi_2_OH3
      integer:: nSi_3_OH,nSi_3_OH2,nSi_3_OH3
      integer:: nSi_4_OH,nSi_4_OH2,nSi_4_OH3
      integer:: nimage,natom_orig,neigh,ncell(125),iOx
      integer,allocatable:: ps_atom(:),del_list(:),ClusterCount(:)
   
      real(wp):: x,y,z,dr(3),nix,niy,niz,boxl_0(3)
      real(wp):: cutoff,cutoff2
      real(wp):: t(0:10)
      real(wp):: ra(3),rb(3),rc(3),r_ba(3),r_ca(3),mag_crossp_3d,temp(3)
      real(wp),allocatable:: ps_rxyz(:,:)
   
      logical:: consistent
      logical,allocatable:: dl_logical(:)
   
      character(len=132):: finput1,finput2,finput3,finput4,carg
      character(len=32):: ctmp,word,atomname,c5
   
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      if(narg /= 2) then
         write(*,*)'usage :'
         write(*,'(/)')
         write(*,*)'  executable_name phase_sep.xyz amorphous.xyz amorphous.conn'
         write(*,'(/)')
      end if
   
      call get_command_argument (1, finput1, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
      call get_command_argument (2, finput2, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
   
      call get_command_argument (3, finput3, leng, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
   
   
      cutoff = 6.0_wp
      cutoff2 = cutoff**2
   
      ! read in the phase separation data
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput1),STATUS='old')
         read(iu,*) ps_natom
         allocate(ps_rxyz(ps_natom,3),ps_atom(ps_natom))
         ps_rxyz = 0.0_wp
   
         read(iu,*) ctmp
         do i = 1,ps_natom
            read(iu,*) ctmp,x,y,z
            ps_atom(i) = name2atom2(ctmp)
            ps_rxyz(i,1:3) = (/ x,y,z /)
         end do
      CLOSE(iu)
   
   
      ! read in coordinates of file amorphous silica system
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput2),STATUS='old')
         read(iu,*) natom
         allocate(rxyz(2*natom,3),proximity(2*natom,4),&
                  proximity_temp(2*natom,4),atom(2*natom))
         natom_orig=natom
         proximity=0

         read(iu,*) boxl(1),boxl(2),boxl(3)
         do i = 1,natom
            read(iu,*) ctmp,x,y,z
            atom(i)=name2atom(ctmp)
            rxyz(i,1:3) = (/ x,y,z /)
         end do
      CLOSE(iu)
   
      call get_free_file_unit(iu)
      OPEN(UNIT=iu,FILE=trim(finput3),STATUS='old')
         call read_conn_pdb(iu)
      CLOSE(iu)
   
      natom_max = 2*natom
   
      call prox_org
      
      
      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*,'pre-crossbond'
      print*, 'not consistent-wrong atoms bonded'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if
      
      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if
   
      boxl_0(1:3)=boxl(1:3)
      boxl2(1:3)=boxl(1:3)/2.0_wp
      natom_max = 2*natom
      allocate(del_list(2*natom))
      del_list=0
      call INIT_NLIST(-boxl2,boxl2,5.0_wp/angstrom,natom_max)
      call new_nlist(1,natom_max)
      call cpu_time(t(0))
      ndel = 0
      del_list = 0
      allocate(dl_logical(natom_max))
      dl_logical(1:natom_max)=.false.
   
      ! loop which marks atoms that overlap the 'boron' phase as .TRUE.
      ! for delete
      atom_loop: do i = 1,ps_natom
         !if (mod(i,100)==0) print*,'i = ', i
         if (ps_atom(i) == iBorate) then
            ic = CELL(ps_rxyz(i,1:3))
            call NEIGCELL(ic,1,neigh,ncell)
            cell_loop: do jj=1,neigh
               nc = ncell(jj)
               if (nc == 0) cycle cell_loop
               j = HOC(nc)
               cell_atom_loop: do while (j /= 0)
                  if(atom(j)==iSilicon) then
                     dr(1:3) = ps_rxyz(i,1:3) - rxyz(j,1:3)
                     call pbc(dr)
                     if ( dot_product(dr,dr) < cutoff2 )  then
                        dl_logical(j)=.TRUE.
                     end if
                  end if
1020              j = LL(j)
               end do cell_atom_loop
            end do cell_loop
         end if
      end do atom_loop
   
   
   
      ! do a simple imaging of the cube in the negative z direction
      ! first move initial centred on origin up a length boxl2(3)
      nimage=2
      rxyz(1:natom,3) = rxyz(1:natom,3) + boxl2(3)
      CALL init_crossbond(nimage,boxl)
      rxyz(natom+1:2*natom,1) = rxyz(1:natom,1)
      rxyz(natom+1:2*natom,2) = rxyz(1:natom,2)
      rxyz(natom+1:2*natom,3) = rxyz(1:natom,3) - boxl(3)
      atom(natom+1:2*natom) = atom(1:natom)
      where (proximity(1:natom,:)/=0) proximity(natom+1:2*natom,:) = &
                      proximity(1:natom,:) + natom
   
      ncross=2*ncross
      natom=2*natom
      nix=1.0_wp
      niy=1.0_wp
      niz=2.0_wp
      CALL crossbond_correct(nix,niy,niz,boxl_0,natom_orig)
   
      call prox_org
   
      ! a check to make sure the imaging and thus connectivity corrections
      ! were done correctly   
      
      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*,'pre-crossbond'
      print*, 'not consistent-wrong atoms bonded'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if
      
      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
      print*, 'not consistent-wrong atoms bonded - post crossbond'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if


      ! delete atoms
      print*,' ndel = ', count(dl_logical(1:natom)==.TRUE.)
      do i = natom,1,-1
         if (dl_logical(i)) call delete_atom(i)
      end do

call cpu_time(t(1))
print*, 'delete time (s): ', t(1)-t(0)

! Find the largest cluster of atoms
      call Init_HKNonLattice(natom)
      atomL = 0
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print '(/a,i8)','n_cluster = ',n_cluster
      call cluster_count(natom,n_cluster,atomL,MaxCluster,ClusterCount)

      print '(/a,i8)','n_cluster = ',n_cluster
      print *,'MaxCluster = ',MaxCluster
      print *,'ClusterCount(MaxCluster) = ',ClusterCount(MaxCluster)
      print *,'sum(ClusterCount) = ',sum(ClusterCount)
      print *,'natom = ',natom

! Delete all atoms not in Largest cluster
      j = natom
      do
         if (atomL(j) /= MaxCluster) then
            atomL(j) = atomL(natom)
            call delete_atom(j)
         end if
         j = j - 1
         if (j == 0) exit
      end do

! convert singly bonded oxygens to OH
      do i = 1,natom
         if (atom(i) == iOxygen) then
            if (count(proximity(i,:) > 0) < 2) atom(i) = iOxygenH
         end if        
      end do
   
! Delete all Si(OH)3 groups
      i = natom
      do
!print '(5i7)',i
!call write_atom_info(6,i)
         if (atom(i) == iSilicon) then
         if (noh(i) == 3) then
            call find_1st_atomtype(i,iOxygen,iOx)
            atom(iOx) = iOxygenH
            call delete_group0(i,iOx)
            i = natom + 1
         end if
         endif
         i = i - 1
         if (i <= 0) exit
      end do

!      ! remove floating atoms
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!   
!     ! sort and delete atoms
!      call shell(floater,del_list)
!      do i = floater,1,-1
!         call delete_atom(del_list(i))
!      end do
!   
!print*, 'floater = ', floater
!print*, 'natom pre-float del',natom
!
!      call re_init_nlist(boxl,5.0_wp)
!      call new_nlist
!
!      !! count number of silicons which are not 4 coordinated
!      !! and which have no OH group
!      ucoord1 = 0
!      ucoord2 = 0
!      ucoord3 = 0
!      coord_loop: do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!            do j=1,4
!               if (proximity(i,j)/=0) then
!                  if (atom(proximity(i,j))==iOxygenH) CYCLE coord_loop
!               end if
!            end do
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               ucoord1 = ucoord1 + 1
!            case(2)
!               ucoord2 = ucoord2 + 1
!            case(3)
!               ucoord3 = ucoord3 + 1
!            end select
!         end if
!         end if
!      end do coord_loop
!      print*, 'ucoord (Si_0_1) = ', ucoord1
!      print*, 'ucoord (Si_0_2) = ', ucoord2
!      print*, 'ucoord (Si_0_3) = ', ucoord3
!      print'(2/)'
!   
!      call prox_org
!
!do
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !                                                      !
!      !              remove any 1 coord Si                   !
!      !              remove any 2 coord Si                   !
!      !           add OH to any 3 coord Si                   !
!      !                                                      !
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      ndel = 0
!      del_list = 0
!      nat_temp = natom
!      at_loop: do i = 1,nat_temp
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!   
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               do j=1,ncmax(atom(i))
!                  if (proximity(i,j)==0) cycle
!                  if (atom(proximity(i,j))==iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!               ndel = ndel + 1
!               del_list(ndel) = i
!   
!            case(2)
!               do j=1,ncmax(atom(i))
!                  if (proximity(i,j)==0) cycle
!                  if (atom(proximity(i,j))==iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!               ndel = ndel + 1
!               del_list(ndel) = i
!   
!            case(3)
!               ! add on an OH to the empty spot
!               natom = natom + 1
!               ra(1:3) = rxyz(proximity(i,1),:)
!               rb(1:3) = rxyz(proximity(i,2),:)
!               rc(1:3) = rxyz(proximity(i,3),:)
!               r_ba = ra - rb
!               call pbc(r_ba)
!               r_ca = ra - rc
!               call pbc(r_ca)
!               mag_crossp_3d = sqrt(dot_product(crossp_3d(r_ba,r_ca),crossp_3d(r_ba,r_ca)))
!               temp = crossp_3d(r_ba,r_ca)/mag_crossp_3d
!               rxyz(natom,:) = rxyz(i,:) + temp*bondl_SiO
!               call pbc(rxyz(natom,:))
!               atom(natom) = iOxygenH
!               proximity(natom,1) = i
!               proximity(i,4) = natom
!            end select
!         end if
!         end if
!      end do at_loop
!   
!      call shell(ndel,del_list)
!   
!print*, 'ndel = ', ndel
!print*, 'natom pre-treat - del',natom
!
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-treat - del',natom
!print'(2/)'
!
!      call prox_org
!   
!      ! remove floating atoms
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!   
!     ! sort and delete atoms
!      call shell(floater,del_list)
!      do i = floater,1,-1
!         call delete_atom(del_list(i))
!      end do
!
!print*, 'floater = ', floater
!
!      CALL Init_HKNonLattice(natom)
!      CALL HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
!      del_list=0
!      ndel=0
!      do i=1,natom
!         if (atomL(i)/=MaxCluster) then
!            ndel=ndel+1
!            del_list(ndel)=i
!         end if
!      end do
!
!      call shell(ndel,del_list)
!      do i=ndel,1,-1
!         call delete_atom(del_list(i))
!      end do
!      deallocate(ClusterCount)
!
!
!      !! count number of silicons which are not 4 coordinated
!      !! and which have no OH group
!      ucoord1 = 0
!      ucoord2 = 0
!      ucoord3 = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) > 0) then
!!            if (ANY(atom(proximity(i,:))==iOxygenH)) CYCLE
!            do j=1,4
!               if (proximity(i,j)/=0) then
!                  if (atom(proximity(i,j))==iOxygenH) CYCLE
!               end if
!            end do
!   
!            select case(count(proximity(i,:)/=0))
!            case(1)
!               ucoord1 = ucoord1 + 1
!            case(2)
!               ucoord2 = ucoord2 + 1
!            case(3)
!               ucoord3 = ucoord3 + 1
!            end select
!         end if
!         end if
!      end do
!print*, 'ucoord (Si_0_1) = ', ucoord1
!print*, 'ucoord (Si_0_2) = ', ucoord2
!print*, 'ucoord (Si_0_3) = ', ucoord3
!print'(2/)'
!
!      if (ucoord1==0 .AND. ucoord2==0 .AND. ucoord3==0) EXIT
!end do
!
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!
!print*, 'floater = ', floater
!print*, 'natom pre-float del',natom
!
!      do k=floater,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-float del',natom
!
!     call prox_org
! 
!     ! go over and delete any clusters of floating atoms
!
!      CALL Init_HKNonLattice(natom)
!      CALL HKNonLattice(natom,proximity,n_cluster,atomL)
!      CALL cluster_count(natom,n_cluster,MaxCluster,ClusterCount)
!      del_list=0
!      ndel=0
!      do i=natom,1-1
!         if (atomL(i)/=MaxCluster) then
!            ndel=ndel+1
!            del_list(ndel)=i
!         end if
!      end do
!
!      do i=1,ndel
!         call delete_atom(del_list(i))
!      end do
!
!do
!      !! delete any FULLY bonded si atoms which have an 3 OH's
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 0) then
!            select case(count(atom(proximity(i,:))==iOxygenH))
!            case(1)
!               ! do nothing
!            case(2)
!               ! do nothing
!            case(3)
!            ndel = ndel + 1
!            del_list(ndel) = i
!   
!            do j = 1,ncmax(atom(i))
!               if (proximity(i,j)==0) stop 'screwed up'
!               if (atom(proximity(i,j)) == iOxygenH) then
!                  ndel = ndel + 1
!                  del_list(ndel) = proximity(i,j)
!               else if (atom(proximity(i,j)) == iOxygen) then
!                  atom(proximity(i,j)) = iOxygenH
!               end if
!            end do
!            case default
!               cycle
!            end select
!         end if
!         end if
!      end do
!   
!      call shell(ndel,del_list)
!   
!      do k=ndel,1,-1
!        call delete_atom(del_list(k))
!      end do
!   
!      !! delete any TRIPLE bonded si atoms which have an 3 OH's
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 1) then
!            ndel = ndel + 1
!            del_list(ndel) = i
!            select case(count(atom(proximity(i,:))==iOxygenH))
!            case(1)
!               ! do nothing
!            case(2)
!               ! no nothing
!            case(3)
!               do j = 1,ncmax(atom(i))
!                  if (proximity(i,j)==0) CYCLE
!                  if (atom(proximity(i,j)) == iOxygenH) then
!                     ndel = ndel + 1
!                     del_list(ndel) = proximity(i,j)
!                  else if (atom(proximity(i,j)) == iOxygen) then
!                     atom(proximity(i,j)) = iOxygenH
!                  end if
!               end do
!            end select
!   
!         end if
!         end if
!      end do
!   
!      call shell(ndel,del_list)
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!print*, 'natom post-triple - del',natom
!      
!      !! convert all singly bonded oxygens to
!      do i = 1,natom
!         if (atom(i)/=iOxygen) CYCLE
!         if (count(proximity(i,:)==0) == 3) then
!            atom(i) = iOxygenH
!         end if
!      end do
!
!print*, 'post convert'
!
!      !! delete any double bonded si atoms which have an OH
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i) == iSilicon) then
!         if (count(proximity(i,:)==0) == 2) then
!            do k=1,4
!               if (proximity(i,k)/=0) then
!                  if (atom(proximity(i,k))==iOxygenH) then
!                     ndel = ndel + 1
!                     del_list(ndel) = i
!                     do j = 1,ncmax(atom(i))
!                        if (proximity(i,j)==0) cycle
!                        if (atom(proximity(i,j)) == iOxygenH) then
!                           ndel = ndel + 1
!                           del_list(ndel) = proximity(i,j)
!                        else if (atom(proximity(i,j)) == iOxygen) then
!                           atom(proximity(i,j)) = iOxygenH
!                        end if
!                     end do
!                  end if
!               end if
!            end do
!         end if
!         end if
!      end do
!
!      call shell(ndel,del_list)
!
!print*, 'ndel = ', ndel
!print*, 'natom double-del',natom
!   
!      do k=ndel,1,-1
!         call delete_atom(del_list(k))
!      end do
!      !! delete any singly bonded si atoms
!      ndel = 0
!      del_list = 0
!      do i = 1,natom
!         if (atom(i)==iSilicon) then
!         if (count(proximity(i,:)==0) == 3) then
!            ndel = ndel + 1
!            del_list(ndel) = i
!            do j = 1,ncmax(atom(i))
!               if (proximity(i,j)==0) cycle
!               if (atom(proximity(i,j)) == iOxygen) then
!                  atom(proximity(i,j)) = iOxygenH
!               end if
!            end do
!         end if
!         end if
!      end do
!
!      !! analyse the bonds in the system
!      nSi_1_OH = 0
!      nSi_1_OH2 = 0
!      nSi_1_OH3 = 0
!      nSi_2_OH = 0
!      nSi_2_OH2 = 0
!      nSi_2_OH3 = 0
!      nSi_3_OH = 0
!      nSi_3_OH2 = 0
!      nSi_3_OH3 = 0
!      nSi_4_OH = 0
!      nSi_4_OH2 = 0
!      nSi_4_OH3 = 0   
!      do i=1,natom
!         if (atom(i) == iSilicon) then
!            n_OH = 0
!            do j = 1,ncmax(atom(i))
!               if (atom(proximity(i,j)) == iOxygenH) then
!                  n_OH = n_OH + 1
!               end if
!            end do
!            coord_type = count(proximity(i,:)/=0)
!            select case(coord_type)
!            case(1)
!               select case(n_OH)
!               case(1)
!                  nSi_1_OH = nSi_1_OH + 1
!               case(2)
!                  nSi_1_OH2 = nSi_1_OH2 + 1
!               case(3)
!                  nSi_1_OH3 = nSi_1_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(2)
!               select case(n_OH)
!               case(1)
!                  nSi_2_OH = nSi_2_OH + 1
!               case(2)
!                  nSi_2_OH2 = nSi_2_OH2 + 1
!               case(3)
!                  nSi_2_OH3 = nSi_2_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(3)
!               select case(n_OH)
!               case(1)
!                  nSi_3_OH = nSi_3_OH + 1
!               case(2)
!                  nSi_3_OH2 = nSi_3_OH2 + 1
!               case(3)
!                  nSi_3_OH3 = nSi_3_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            case(4)
!               select case(n_OH)
!               case(1)
!                  nSi_4_OH = nSi_4_OH + 1
!               case(2)
!                  nSi_4_OH2 = nSi_4_OH2 + 1
!               case(3)
!                  nSi_4_OH3 = nSi_4_OH3 + 1
!!               case(4)
!!                  print*, 'ERROR: floating silicon'
!               end select
!            end select
!         else if (atom(i) == iOxygen) then
!            ! do nothing
!         else if (atom(i) == iOxygenH) then
!            ! do nothing
!         else
!            print*, 'error in atom types'
!            print*, 'atom(',i,') = ', atom_name(atom(i))
!            print'(2/)'
!            print*, 'program terminated prematurely'
!            print'(/)'
!            stop
!         end if
!      end do
!   
!      !! bonding summary
!      print'(2/)'
!!      print*, 'nSi_0H   = ', nSi_OH
!!      print*, 'nSi_0H2  = ', nSi_OH2
!!      print*, 'nSi_0H3  = ', nSi_OH3
!      print*, 'Singly bonded silicons'
!      print*, 'nSi_1_OH  = ',  nSi_1_OH
!      print*, 'nSi_1_OH2 = ',  nSi_1_OH2
!      print*, 'nSi_1_OH3 = ',  nSi_1_OH3
!      
!      print'(2/)'
!      print*, 'Double bonded silicons'
!      print*, 'nSi_2_OH  = ',  nSi_2_OH
!      print*, 'nSi_2_OH2 = ',  nSi_2_OH2
!      print*, 'nSi_2_OH3 = ',  nSi_2_OH3
!      
!      print'(2/)'
!      print*, 'Triple bonded silicons'
!      print*, 'nSi_3_OH  = ',  nSi_3_OH
!      print*, 'nSi_3_OH2 = ',  nSi_3_OH2
!      print*, 'nSi_3_OH3 = ',  nSi_3_OH3
!      
!      print'(2/)'
!      print*, 'Fully bonded silicons'
!      print*, 'nSi_4_OH  = ',  nSi_4_OH
!      print*, 'nSi_4_OH2 = ',  nSi_4_OH2
!      print*, 'nSi_4_OH3 = ',  nSi_4_OH3
!      print'(2/)'
!
!      ! if all ducks are in a row break out
!      if (nSi_1_OH==0 .AND.&
!          nSi_1_OH2==0 .AND.&
!          nSi_1_OH3==0 .AND.&
!          nSi_2_OH==0 .AND.&
!          nSi_2_OH2==0 .AND.&
!          nSi_2_OH3==0 .AND.&
!          nSi_3_OH==0 .AND.&
!          nSi_3_OH2==0 .AND.&
!          nSi_3_OH3==0 .AND.&
!          nSi_4_OH3==0) exit
!end do
!   
!   
!      !! convert all singly bonded oxygens to
!      do i = 1,natom
!         if (atom(i)/=iOxygen) CYCLE
!         if (count(proximity(i,:)==0) == 3) then
!            atom(i) = iOxygenH
!         end if
!      end do
!   
!      del_list=0
!      floater = 0
!      do k = 1, natom
!         if (count(proximity(k,:)/=0) == 0) then
!            floater = floater + 1
!            del_list(floater)=k
!         end if
!      end do
!      call shell(floater,del_list)
!   
!   print*, 'floater = ', floater
!   print*, 'natom pre-float del',natom
!   
!      do k=floater,1,-1
!         call delete_atom(del_list(k))
!      end do
   
    
      !! analyse the bonds in the system
      nSi_1_OH = 0
      nSi_1_OH2 = 0
      nSi_1_OH3 = 0
      nSi_2_OH = 0
      nSi_2_OH2 = 0
      nSi_2_OH3 = 0
      nSi_3_OH = 0
      nSi_3_OH2 = 0
      nSi_3_OH3 = 0
      nSi_4_OH = 0
      nSi_4_OH2 = 0
      nSi_4_OH3 = 0
      do i=1,natom
         if (atom(i) == iSilicon) then
            n_OH = 0
            do j = 1,ncmax(atom(i))
               if (atom(proximity(i,j)) == iOxygenH) then
                  n_OH = n_OH + 1
               end if
            end do
      
      
            coord_type = count(proximity(i,:)/=0)
            select case(coord_type)
            case(1)
               select case(n_OH)
               case(1)
                  nSi_1_OH = nSi_1_OH + 1
               case(2)
                  nSi_1_OH2 = nSi_1_OH2 + 1
               case(3)
                  nSi_1_OH3 = nSi_1_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(2)
               select case(n_OH)
               case(1)
                  nSi_2_OH = nSi_2_OH + 1
               case(2)
                  nSi_2_OH2 = nSi_2_OH2 + 1
               case(3)
                  nSi_2_OH3 = nSi_2_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(3)
               select case(n_OH)
               case(1)
                  nSi_3_OH = nSi_3_OH + 1
               case(2)
                  nSi_3_OH2 = nSi_3_OH2 + 1
               case(3)
                  nSi_3_OH3 = nSi_3_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            case(4)
               select case(n_OH)
               case(1)
                  nSi_4_OH = nSi_4_OH + 1
               case(2)
                  nSi_4_OH2 = nSi_4_OH2 + 1
               case(3)
                  nSi_4_OH3 = nSi_4_OH3 + 1
!               case(4)
!                  print*, 'ERROR: floating silicon'
               end select
            end select
      
      
         else if (atom(i) == iOxygen) then
            ! do nothing
      
         else if (atom(i) == iOxygenH) then
            ! do nothing
      
         else
            print*, 'error in atom types'
            print*, 'atom(',i,') = ', atom_name(atom(i))
            print'(2/)'
            print*, 'program terminated prematurely'
            print'(/)'
            stop
         end if
      end do
         
      !! bonding summary
      print'(2/)'
      ! print*, 'nSi_0H   = ', nSi_OH
      ! print*, 'nSi_0H2  = ', nSi_OH2
      ! print*, 'nSi_0H3  = ', nSi_OH3
      print*, 'Singly bonded silicons'
      print*, 'nSi_1_OH  = ',  nSi_1_OH
      print*, 'nSi_1_OH2 = ',  nSi_1_OH2
      print*, 'nSi_1_OH3 = ',  nSi_1_OH3
      
      print'(2/)'
      print*, 'Double bonded silicons'
      print*, 'nSi_2_OH  = ',  nSi_2_OH
      print*, 'nSi_2_OH2 = ',  nSi_2_OH2
      print*, 'nSi_2_OH3 = ',  nSi_2_OH3
      
      print'(2/)'
      print*, 'Triple bonded silicons'
      print*, 'nSi_3_OH  = ',  nSi_3_OH
      print*, 'nSi_3_OH2 = ',  nSi_3_OH2
      print*, 'nSi_3_OH3 = ',  nSi_3_OH3
      
      print'(2/)'
      print*, 'Fully bonded silicons'
      print*, 'nSi_4_OH  = ',  nSi_4_OH
      print*, 'nSi_4_OH2 = ',  nSi_4_OH2
      print*, 'nSi_4_OH3 = ',  nSi_4_OH3
      print'(2/)'
      
      call prox_org
      
!      do i=1,natom
!      if (count(proximity(i,:)/=0)<ncmax(atom(i))) then
!      if (atom(i)==iOxygenH) then
!      if (count(proximity(i,:)/=0)==1) cycle
!      end if
!            print*, 'atom(',i,') = ', atom_name(atom(i))
!            print'(a5,1x,i6,a6,1x,4i6)', 'prox(',i,',:) = ',(proximity(i,j),j = 1,ncmax(atom(i)))
!            print'(2/)'
!         end if
!      end do
!
      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
         print*, 'not consistent 2 - end'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if
      
      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
      print'(2/)'
         print*, 'not consistent 1 - end'
      print'(/)'
      call write_atom_info(6,ifirst)
      print'(3/)'
      stop
      end if

   call get_free_file_unit(iu)
   OPEN(UNIT=iu,FILE='large_etch.xyz',STATUS='UNKNOWN')
      write(iu,*) natom
      write(iu,*) 'Etched_Structure'
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)), rxyz(i,1:3)
      end do
   CLOSE(iu)

   call get_free_file_unit(iu)
   OPEN(UNIT=iu,FILE='large_etch.conn',STATUS='UNKNOWN')
      write(iu,*) natom
      write(iu,*) 'connectivity data for amorphous.xyz'
      do i = 1,natom
         write(iu,'(5i7)') i,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
   CLOSE(iu)

110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3)



CONTAINS

   SUBROUTINE prox_org
      ! arranges the non-zero elements of each proximity
      ! array from the leftmost side
      USE connectivity_mod, only: proximity
      USE global_vars_mod, only: natom
      IMPLICIT NONE
      integer:: i,ptmp(4),nat
      do i=1,natom
         ptmp(1:4) = proximity(i,1:4)
         proximity(i,:) = 0
         nat = 0
         do j=1,4
            if (ptmp(j)/=0) then
               nat = nat + 1
               proximity(i,nat) = ptmp(j)
            end if
         end do
      end do
   END SUBROUTINE prox_org



   pure function name2atom2(c)
      integer:: name2atom2
      character(*),intent(in):: c
      integer:: i
      do i = 0,1
         if(trim(c) == trim(atom_name2(i)))then
            name2atom2 = i
            exit
         end if
      end do
      if (i > ntyp) name2atom2 = -1
   end function

   SUBROUTINE read_conn_pdb(iu)
      USE connectivity_mod
      USE atom_types_mod
      integer,intent(in):: iu
      integer:: nat,i,j,iatom
      character(len=132):: line
      character(len=80):: ctmp
      character(len=32):: c5

      ! read in a formatted connectivity file
      read(iu,*) nat
      read(iu,*) ctmp
      do i=1,nat
         line=''
         read(unit=iu,fmt='(a132)',iostat=ios) line
         c5=line(1:6)
         read(unit=c5,fmt=*)iatom
         do j=1,ncmax(atom(iatom))
            c5=line(6*(j)+1:6*(j)+6)
            read(unit=c5,fmt=*) proximity(iatom,j)
         end do
      end do

   END SUBROUTINE read_conn_pdb


END PROGRAM etch_amorphous

