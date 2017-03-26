include 'precision_mod.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'bond_list_mod.f90'
include 'crossbond_mod_v4.f90'
include 'frames_mod.f90'

PROGRAM crossbond_test

! the aim of this program is to test the crossbond module that 
! has been written. the simbox_mod will be used to create an
! initial box using the silica unit cell. the conectivity etc. 
! will then be calculated using the brute force method. the 
! neighbour list method being used currently will be neglected 
! for testing purposes. 
! using the newly imaged box it will be imaged again the same 
! number of times to give a 3x3x3 large cube of ~ 17500 atoms
! the connectivity and atom numbers should also be imaged in 
! this process. using this process the only connectivity issue
! will be at the boundaries and this is the point where the new 
! code will be tested. 

! ACT I - box creation
      USE precision_mod, only: wp      
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE frames_mod
      USE atom_list_mod
      USE crossbond_mod
      USE files_mod
      
      IMPLICIT NONE
      integer,parameter:: natom_unit = 24
      integer,parameter:: nl = 1 
      integer,parameter:: nunitc = 2*nl + 1
      real(wp):: r3(3),rj(3),bondl_orig=1.55_wp,bondl_scale
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)!,ia(:)
      real(wp):: CL,x,y,z,dr(3),rverlet!,bl1,bl2
      integer:: j,i,n,ib,k,imve,ix,iy,iz,iu3,ic,jj,nc!,natom
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_crossbondtot
      character(32):: ctmp
!      type(atom_list):: lst_positive_x,lst_negative_x,lstSi
!      type(atom_list):: lst_positive_y,lst_negative_y
!      type(atom_list):: lst_positive_z,lst_negative_z
      n_crossbondtot = 1000 ! 1000 should cover the num of crossbonds in any dir.
!
!      allocate(lst_negative_x%ind(n_crossbondtot),lst_negative_y%ind(n_crossbondtot))
!      allocate(lst_negative_z%ind(n_crossbondtot))
!      allocate(lst_positive_x%ind(n_crossbondtot),lst_positive_y%ind(n_crossbondtot))
!      allocate(lst_positive_z%ind(n_crossbondtot))
      natom_max = natom_unit*nunitc*nunitc*nunitc
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      nseed = 0
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:2*natom_max))

      OPEN(UNIT=13,FILE='unitcellreal_3.xyz')
      !open(13,file='1.xyz')
      read(13,*) natomc
      read(13,*) ctmp!,ctmp,ctmp,it
      do i = 1,natomc
         read(13,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/angstrom
         atom(i) = name2atom(ctmp)
      end do
      CLOSE(13,STATUS='keep')

      
      allocate(rxyzc(1:natomc,3))
      allocate(proximity(natom_max*(2*nl+1)**3,4))
      print*, 'size(proximity(:,1)) = ',size(proximity(:,1))
      allocate(atomc(1:natomc))
!      allocate(ia(10))
      proximity = 0.0_wp
      nattached = 0
      CL = 7.17_wp/angstrom
      boxl = nunitc*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl

      bondl_scale = bondl_SiO/bondl_orig
      CL = CL*bondl_scale
      boxl=CL*nunitc
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl


natom = natomc

! initial connectivity for unit cell
      do k=1,natom
         n=0
         do j=1,natom
            if(j == k) cycle           
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.2_wp*bondl_SiO)**2) then
                n=n+1
                proximity(k,n)=j
               if (atom(k) == atom(j)) stop 'error: 2 like atoms bonded'
            end if
         end do

      end do

      CALL bond_list 

! creating a simulation box
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
     
! write reduced cell contents to a file
      CALL get_free_file_unit(iu3)
      OPEN(UNIT=iu3,FILE='unitcellreal_4.xyz')
      write(iu3,*) natom
      write(iu3,'(a,i6)') 'Silica'
      do i = 1,natom  
         write(iu3,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*angstrom
      end do
      write(iu3,*)
      CLOSE(iu3, status='keep')
      deallocate(nbond,ibond)
      allocate(nbond((ib+1)*natom),ibond(2,natom*(ib+1)*2))

      proximity = 0
      
! connectivity for imaged cell
      do k=1,natom
         n=0
         do j=1,natom
            if(j == k) cycle           
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.2_wp*bondl_SiO)**2) then
                n=n+1
                proximity(k,n)=j
                if (atom(k) == atom(j)) stop 'error: 2 like atoms bonded'
            end if
         end do
!print*, 'proximity(k,:) = ', proximity(k,:)
      end do
 

do k = 1, natom
if (proximity(k,2)==0) then
print*, 'atom(',k,') = ', atom(k)
print*, 'proximity(',k,',:) = ', proximity(k,:)
print*, ''
end if
end do

 
      CALL bond_list
      ! Make arrays of crossing bonds
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i=1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (ABS(r3(1)) > boxl2) then
            n_crossbond_x = n_crossbond_x+1
            crossbond_x(n_crossbond_x)=i
         end if
         if (ABS(r3(2)) > boxl2) then
            n_crossbond_y = n_crossbond_y+1
            crossbond_y(n_crossbond_y)=i
         end if
         if (ABS(r3(3)) > boxl2) then
            n_crossbond_z = n_crossbond_z+1
            crossbond_z(n_crossbond_z)=i
         end if
      end do

      do i = 1,n_crossbond_x
         ib = crossbond_x(i)
         if (rxyz(ibond(1,ib),1) > 0.0) then
            call add_to_list(lst_positive_x,ibond(1,ib))
            call add_to_list(lst_negative_x,ibond(2,ib))
            if(rxyz(ibond(2,ib),1) > 0.0) then
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
            if(rxyz(ibond(2,ib),2) > 0.0) then
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
            if(rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary'
            end if
         else
            call add_to_list(lst_negative_z,ibond(1,ib))
            call add_to_list(lst_positive_z,ibond(2,ib))
         end if
      end do
!
!print*, 'crossbond_x(:) = ', crossbond_x
!print*, 'crossbond_y(:) = ', crossbond_y
!print*, 'crossbond_z(:) = ', crossbond_z
!
! creating a simulation box

print*, 'natom_unit = ', natom_unit

      
      CALL crossbond_correct(nl,natom_unit)
print*, 'hello #############'
      call new_frame_file(imve,'frame',3)
      call write_frame(imve,1,natom)
      
END PROGRAM crossbond_test 
