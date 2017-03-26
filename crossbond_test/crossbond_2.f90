! READ IN THE INITIAL DATA
! WHEN HAVE THE DATA THEN WE NEED TO IDENTIFY ALL ATOMS INVOLVED IN CROSSBONDS
! THIS WILL BE DONE WITH THE CROSSBOND CODE IN DEVELOPEMENT. IT IS NOW NOT
! REQUIRED TO ACTUALLY SEGREGATE THE BONDS INTO THE DIFFERENT DIRECTIONS AS THE
! PROXIMITY VALUES OF EACH WILL NOW BE SET TO ZERO AND CALCULATED USING THE
! NEIGHBOUR LIST ALGORITHM. 
! ALTERNATIVELY TO THE ORIGINAL APPROACH TO THE PROBLEM OF CORRECTING THE
! CROSSBONDS TO THEIR IMAGES IN THE VARIOUS DIRECTIONS, THIS WILL IMAGE ALL ATOM
! COORDINATES AND ALL RESPECTIVE PROXIMITY VALUES EXCEPT THOSE OF ATOMS INVOLVED
! IN CROSSBONDS. IN THE LATTER CASE THE PROXIMITY VALUES OF THE ATOM ARE SET TO
! ZERO. 
! FOLLOWING THIS THE SITUATION EXISTS THAT ALL THE CONNECTING ATOMS HAVE NO
! CONNECTIVITY. THESE ARE LOOPED OVER AND AN ATTEMPT IS MADE TO CORRECT THE
! CONNECTIVITY OF THE SYSTEM. THIS IS DONE USING THE INITIAL METHOD OF ATOMS
! BEING WITHIN AN SI-O BOND LENGTH AWAY FROM EACH OTHER. IF IT HAPPENS THAT ONE
! OF THESE ATOMS HAS A NUMBER OF NEIGHBOURS CLOSER THAN THE CUT-OFF DISTANCE
! WHICH IS HIGHER THAN THE DESIRED COORDINATION THEN WE GO BACK TO THE SAME ATOM
! IN ORIGINAL BOX AND TAKE THE "EXACT" DISTANCES BETWEEN IT AND IT'S NEIGHBOURS.
! USING THESE DISTANCES WE SEARCH FOR NEIGHBOURS WITH CORRESPONDING EXACT
! DISTANCES FROM TARGET ATOM. 

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



! READ IN THE INITIAL DATA OF STRUCTURE
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
   IMPLICIT NONE

   ! variable declaration
   integer,parameter:: nl = 1 
   integer,parameter:: nunitc(3) = 2*nl + 1
   integer:: ip,itmp,ifirst
   integer:: i,ii,j,jj,k,kk,m,n !counters and dummy args
   integer:: iu,ios,ia,ig,ilast,iatom,natom_orig,ix,iy,iz,ib
   integer:: n_crossbond,ncross,nc,ic,neigh2,neighmx=125
   integer,allocatable:: atomc(:),proximity_temp(:,:)
   integer,allocatable:: cross_atom(:),atom_temp(:),crossbond(:)
   
   real(wp):: timep(0:10),bond_lens(4)
   real(wp):: r3(3),rj(3),dr(3),bondl_orig=1.55_wp
   real(wp),allocatable:: rxyz_temp(:,:)
   real(wp):: CL,x,y,z,boxl_0(3)

   character(len=32):: carg,word,ctmp
   character(6):: c6
   character(5):: c5

   logical:: consistent

   !detail on the input box's size
   natom = 3000
   natom_max = 8*natom

   CALL cpu_time(timep(0)) 
   
   CALL get_free_file_unit(iu)
   allocate(atom_temp(20000))
   allocate(atom(natom_max))
   allocate(rxyz_temp(natom,3),proximity_temp(natom,4))
   rxyz_temp=0.0_wp
   proximity_temp=0
   
   OPEN(UNIT=iu, FILE='3k_SiO2.pdb',STATUS='old') !change the file name to suit
   read(iu,*) ctmp,CL
   do i=1,natom
      read (iu,*) c6,itmp,ctmp,itmp,(rxyz_temp(i,j),j=1,3)
      atom_temp(i) = name2atom(trim(ctmp))
   end do
   do i=1,natom
      read(iu,'(a32)')ctmp
      do j = 1,ncmax(atom_temp(i))
         c5 = ctmp(6+5*(j)+1:6+5*(j)+5)
         !print*, c5
         read( unit=c5,fmt=* ) proximity_temp(i,j)
      end do 
   end do
   CLOSE(iu, STATUS='KEEP')
!detail on the input box's size
   !CL=CL*(bondl_SiO/bondl_orig)
   boxl(:) = CL
   boxl_0(:) = boxl(:)
   boxl2 = boxl/2.0_wp

   
   natom_orig = natom
   print*, 'natom',natom
   print*, 'natom_max',natom_max
   allocate(rxyz(natom_max,3),proximity(natom_max,4))
   proximity = 0
   proximity(1:natom,:) = proximity_temp(1:natom,:)
   rxyz(1:natom,:) = rxyz_temp(1:natom,:)
   atom(1:natom) = atom_temp(1:natom)

   CALL bond_list
   print*, nbondtot
   
!  CHECK FOR CROSSBONDS. 
   allocate(crossbond(nbondtot))
   n_crossbond = 0
   do i=1,nbondtot ! nbondtot calculated in bond_list
      r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:) 
      r3 = SQRT(dot_product(r3,r3))      
      if (r3(1) > boxl2(1) .OR. r3(2) > boxl2(2) .OR. &
          r3(3) > boxl2(3)) then
         n_crossbond = n_crossbond + 1 
         crossbond(n_crossbond) = i
      end if 
   end do

! store all of the atom numbers that are involved in crossbonds and the
! connectivity of these atoms
   allocate(cross_atom(12*nbondtot*2)) ! 12 images and 2 atoms per bond
   j=1
   ncross = 0
   do i=1,n_crossbond
      cross_atom(j) = ibond(1,crossbond(i))
      cross_atom(j+1) = ibond(2,crossbond(i))
      ncross = ncross + 2
      j=j+2
   end do

   do i =1, ncross
      proximity(cross_atom(i),:) = 0
   end do

   ! initial translation of box occurs as follows 
   rxyz(1:natom,1) = rxyz(1:natom,1) + boxl2(1)
   rxyz(1:natom,2) = rxyz(1:natom,2) + boxl2(2)
   rxyz(1:natom,3) = rxyz(1:natom,3) + boxl2(3)
   
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!change the looping number of layers to give correct structure
!as well as this, may have to 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   ib = 0    
   do iz = -nl,0
   do iy = -nl,0
   do ix = -nl,0
      if((ix==0).and.(iy==0).and.(iz==0)) cycle 
      ib = ib + 1         
! image the coordinates
      rxyz(ib*natom+1:(ib+1)*natom,1) = rxyz(1:natom,1) + boxl(1)*ix
      rxyz(ib*natom+1:(ib+1)*natom,2) = rxyz(1:natom,2) + boxl(2)*iy
      rxyz(ib*natom+1:(ib+1)*natom,3) = rxyz(1:natom,3) + boxl(3)*iz
      
! image the atom types      
      atom(ib*natom+1:(ib+1)*natom) = atom(1:natom)         
      
!! image the proximity array also
      where (proximity(1:natom,:)/=0) proximity(ib*natom+1:(ib+1)*natom,:) = &
                      proximity(1:natom,:) + ib*natom
   end do
   end do
   end do
   ncross = ib*ncross
   natom = natom*(ib+1)
   
! need to specify the changes in the box length directions for initialisation of
! the neighbour list. 
   boxl = (/boxl_0(1)*2.0_wp,boxl_0(2)*2.0_wp,boxl_0(3)*2.0_wp/) 
   boxl2 = boxl/2.0_wp

   call INIT_NLIST(boxl,5.0_wp)
   call NEW_NLIST

! the following is to get the connectivity of the edge atoms using the 
! adjusted nlist
main_loop: do i = 1, natom
      if ( .NOT. ALL( proximity(i,:)==0 ) ) cycle
      r3(1:3) = rxyz(i,1:3)
      ic = CELL(r3)
      CALL NEIGCELL(ic,1,neigh,ncell)
      n = 0
      cell_loop: do jj = 1, neigh
         nc = ncell(jj)
         if(nc==0) cycle cell_loop
         j=HOC(nc)
         cell_atom_loop: do while (j/=0)
            dr=r3-rxyz(j,1:3)
            if (j==i) GOTO 10
            CALL pbc(dr)
            if (dot_product(dr,dr)<(1.1*bondl_SiO)**2) then 
               n=n+1
               proximity(i,n) = j
               if (n > 4) then   
                  print*, 'defect zone,n = ', n
                  print*, 'prox(',i,',:) = ', proximity(i,:)
                  n=0
                  proximity(i,:) = 0
                  if (i>natom_orig) then
                     iatom = mod(i,natom_orig)
                  else 
                     iatom = i
                  end if
                  ! need to have an array containing the connectivity
                  ! of the target atom i. iatom is the original of
                  ! the image under investigation. 
                  ! change boxl to proper value
                  boxl = (/boxl_0(1),boxl_0(2),boxl_0(3)/) 
                  boxl2 = boxl/2.0_wp 
                  
                  bond_lens = 0.0_wp
                  do k = 1, ncmax(atom(iatom))
                     dr = rxyz(iatom,:) - rxyz(proximity_temp(iatom,k),:)
                     CALL pbc(dr)
                     bond_lens(k) = dot_product(dr,dr)
                  end do
                  
                  !reset boxl to proper value
                  boxl = (/boxl_0(1)*2.0_wp,boxl_0(2)*2.0_wp,boxl_0(3)*2.0_wp/) 
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
                           if (ABS( dot_product(dr,dr)-bond_lens(ii) )<=1.0E-8_wp) then
                              EXIT
                           else
                              GOTO 20
                           end if 
                        end do       
                        n=n+1
                        proximity(i,n) = m          
                        if (n>4) stop 'need better condition'
                        if (n==ncmax(atom(i))) cycle main_loop
20                      CONTINUE
                        m=LL(m)
                     end do inner_cell_atom_loop
                  end do inner_cell_loop
                  cycle cell_loop
               end if
            end if
10          CONTINUE
            j=LL(j)
         end do cell_atom_loop
      end do cell_loop 
   end do main_loop
   
   boxl = (/boxl_0(1)*2.0_wp,boxl_0(2)*2.0_wp,boxl_0(3)*2.0_wp/) 
   boxl2 = boxl/2.0_wp  
   call bond_list
   CALL cpu_time(timep(1))
   print*, 'time taken = ', timep(1) - timep(0)
   
   call check_proximity(consistent,ifirst,k)
   if (.NOT. consistent) then
      print*, 'not consistent'
      print*, 'proximity(',ifirst,',:)=', proximity(ifirst,:)
      print*, 'proximity(',k,',:)=', proximity(k,:)
   end if
   
   CALL get_free_file_unit(iu)
   CALL new_frame_file(iu,'imaged',0)
   CALL write_frame(iu,1,natom)
   close(iu)
   
END PROGRAM crossbond_correct

