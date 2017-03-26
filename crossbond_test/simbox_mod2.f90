MODULE simbox_mod
   USE precision_mod, only: wp
   implicit none

CONTAINS
   
   SUBROUTINE SIMBOX(boxl0,rverlet)
      !include 'precision_mod.f90'
      !include 'files_mod.f90'
      !include 'rand_mod.f90'
      !include 'sort_mod.f90'
      !include 'global_vars.f90'
      !include 'seaton_mod.f90'
      !include 'constants_mod.f90'
      !include 'atom_types.f90'
      !include 'coordinates_mod.f90'
      !include 'connectivity_mod.f90'
      !include 'bond_list_mod.f90'
      !include 'frames_mod.f90'

      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE frames_mod
      USE nlist_mod
      USE verlet_list
      implicit none
      real(wp),intent(in):: rverlet,boxl0
      integer,parameter:: natom_unit = 24
      integer,parameter:: nl = 3 
      integer,parameter:: nunitc = 2*nl + 1
      real(wp):: r3(3),rj(3)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:),ia(:)
      real(wp):: CL,x,y,z!,bl1,bl2
      integer:: j,i,n,ib,k,imve,natomc,ix,iy,iz,iu3
      integer:: red_cell
      character(32):: ctmp

      open(unit=18,file='neighbour.out')

      natom_max = natom_unit*nunitc*nunitc*nunitc
      natomc = natom_unit
!      print*, 'natom_max = ', natom_max 
!      print*, 'natomc = ', natomc
!      print*, 'nunitc = ', nunitc
      nseed = 0
      allocate(rxyz(1:natom_max,3))
      allocate(rxyzc(1:natom_unit,3))
      allocate(atom(1:natom_max))
      allocate(proximity(natom_max,4))
      allocate(atomc(1:natom_unit))
      allocate(ia(10))
      proximity = 0
      nattached = 0
      boxl2 = boxl0/2.0_wp
      boxli = 1.0_wp/boxl0
      CL = 7.17_wp/angstrom
      
!print*,'rxyz(1,:) = ',size(rxyz(:,1))
!print*,'rxyzc(1,:) = ',size(rxyzc(:,1))
!print*,'atom(:) = ',size(atom(:))
!print*,'proximity(1,:) = ',size(proximity(:,1))
!print*,'atomc(:) = ',size(atomc(:))

! Initialise the neighbour list and verlet list
      CALL INIT_NLIST(boxl0,boxl0,boxl0,5.0_wp/angstrom)

      CALL Init_Verlet_List(rverlet,2.6_wp)
      
      OPEN(UNIT=13,FILE='unitcellreal.xyz')
      !open(13,file='1.xyz')
      read(13,*) natomc
      read(13,*) ctmp!,ctmp,ctmp,it
      do i = 1,natomc
         read(13,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/angstrom
!         print*, 'atom type = ', ctmp
!         print*, 'rxyz(',i,',:) = ',rxyz(i,1:3)
         
         atom(i) = name2atom(ctmp)
      end do
      CLOSE(13,STATUS='keep')
           
! initialise the connectivity for unit cell
      do k=1,natomc
         n=0
         do j=1,natomc
            if(j == k) cycle           
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.1_wp*bondl_SiO)**2) then
                n=n+1
                proximity(k,n)=j
            end if
         end do
      end do

      rxyzc(1:natom_unit,:) = rxyz(1:natom_unit,:)
      atomc(1:natom_unit) = atom(1:natom_unit)
      rj(:) = 0.5_wp*(rxyz(11,:)+rxyz(18,:))
      do i =1, size(rxyz(:,1))
         rxyz(i,:) = rxyz(i,:) - rj(:)
      end do

! write translated contents to a file
      OPEN(UNIT=16,FILE='unitcellreal_2.xyz')
      write(16,*) natomc   
      write(16,'(a,i6)') 'Silica'
      do i = 1,natomc  
         write(16,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*angstrom
      end do
      write(16,*)
      CLOSE(16, status='keep')

!create a temporary cell that is used to store the new atom coords after
!translation and deletion of atoms
! ia are the atoms which i want removed from cell
      ia = (/3,7,8,15,16,18,25,27,31,34/)!
      natom = natomc
      do i =size(ia),1,-1
!         print*, 'atom_type', atom_name(atom(ia(i)))
         call delete_atom(ia(i))
      end do

!-----------------------------------------------------------------------------
! create a new neigh_list
      CALL NEW_NLIST
!-----------------------------------------------------------------------------
      
! write reduced cell contents to a file
      OPEN(UNIT=14,FILE='unitcellreal_3.xyz')
      red_cell = natomc-size(ia)
      write(14,*) red_cell
      write(14,'(a,i6)') 'Silica'
      do i = 1,red_cell  
         write(14,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*angstrom
      end do
      write(14,*)
      CLOSE(14, status='keep')
!------------------------------------------------------------

!      boxl = nunitc*CL
!print*, ' boxl =', boxl
!      boxl2 = boxl/2.0_wp
!      boxli = 1.0_wp/boxl

! creating a simulation box
natomc = natom_unit
!natomc = 1
      ib = 0
      do iz = -nl,nl
      do iy = -nl,nl
      do ix = -nl,nl
         if((ix==0).and.(iy==0).and.(iz==0))cycle
         ib = ib + 1
         rxyz(ib*natomc+1:(ib+1)*natomc,1) = rxyz(1:natomc,1) + CL*ix
         rxyz(ib*natomc+1:(ib+1)*natomc,2) = rxyz(1:natomc,2) + CL*iy
         rxyz(ib*natomc+1:(ib+1)*natomc,3) = rxyz(1:natomc,3) + CL*iz
         atom(ib*natomc+1:(ib+1)*natomc) = atom(1:natomc)
      end do
      end do
      end do
      natom = natom_unit*nunitc*nunitc*nunitc


      CALL NEW_NLIST     
      
! the following is to get the connectivity using the nlist
   
      do i = 1, natom
         r3(1:3) = rxyz(i,1:3)
         ic = CELL(r3)
         call NEIGCELL(ic,1,neigh,ncell)
      cell_loop: do jj = 1, neigh
            n=0
            nc = ncell(jj)
            if(nc==0) cycle cell_loop
            j=HOC(nc)
            cell_atom_loop: do while (j/=0)
               dr=r3-rxyz(j,1:3)
               if (dr==0) then cycle cell_loop
               call pbc(dr)
               if(dot_product(dr,dr)<(1.1*bondl_SiO)**2) then 
                  n=n+1
                  proximity(i,n) = j
if (n>4) then
print*, 'n = ', n
print*, 'prox = ', proximity(i,:)
end if
               end if
               j=LL(j)
            end do 
         end do cell_loop
      end do 

      
!      ! determine connectivity
!      j=1
!      do k=1,natom
!         n=0
!         do j=1,natom
!            if(j == k) cycle
!            if(atom(k) == atom(j)) cycle
!            r3(:)=rxyz(j,:)-rxyz(k,:)
!            CALL pbc(r3)
!            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
!            if (dot_product(r3,r3) < (1.1_wp*bondl_SiO)**2) then
!                n=n+1
!                proximity(k,n)=j
!            end if
!         end do
!         
!!if (mod(k,1000)==0) then
!!print*, 'k =', k
!!end if
!
!      end do
  
      do i=1,natom
         write (18,'(a2,4(1x,i14.8))') atom_name(atom(i)),proximity(i,1:4)
      end do
      close(18,status='KEEP')

      ! write complete xyz data of cell to a file
      open(UNIT= 15, FILE='cell.xyz')
      call write_xyz(15,0)
      close(15,STATUS='keep')
      
      open(unit=19,file='prox_data.out')   
!do i=1,natom
!print*, 'proximity(',i,',:)', proximity(i,:)
!end do

      call write_proximity(19,natom)
      CLOSE(19,STATUS='keep')
      
      call new_frame_file(imve,'frame',0)
      call write_frame(imve,0,natom)
     
   END SUBROUTINE SIMBOX

   
END MODULE simbox_mod


