
include 'precision_mod.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'HKNonLattice2.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
!include 'atom_list_mod.f90'
include 'bond_list_mod.f90'
include 'frames_mod.f90'


PROGRAM SIMBOX
      USE precision_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
!      USE atom_list_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE HKNonLattice_mod
      USE frames_mod
      implicit none
      integer,parameter:: natom_unit = 24
      integer,parameter:: nunitc = 1
      real(wp):: r3(3)
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp):: CL,bl1,bl2,r3l
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic=0,natomc,ix,iy,iz
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
!      type(atom_list):: lst_positive_x,lst_negative_x,lstSi
!      type(atom_list):: lst_positive_y,lst_negative_y
!      type(atom_list):: lst_positive_z,lst_negative_z
!
!      open(unit=16,file='connect_struct.out')
!      open(unit=14,file='simbox.out')
      open(unit=15,file='neighbour.out')
      natom_max = natom_unit*nunitc*nunitc*nunitc
      natomc = natom_unit
      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(proximity(natom_max,4))
      proximity = 0
      print*, 'natom_max = ', natom_max
      print*, 'natom_unit = ', natom_unit


! read input data
      !should change the input filename to the realistic data if that is
      !required
      open(unit=13,FILE='unitcellreal_1.dat')
      do i=1,natom_unit
         read(13,*) rxyz(i,1),rxyz(i,2),rxyz(i,3)
         if(i <= 8) then
            atom(i) = iSilicon
         else
            atom(i) = iOxygen
         end if
      end do
      close(13,status='keep')

!------------------------------------------------------------
      CL = 7.17_wp
      boxl = nunitc*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!      bl1 = CL*sqrt(2.0_wp)/8.0_wp
!      bl2 = CL*sqrt(1.5_wp)/4.0_wp
!
! Creating the simulation box
!
!      natom = 0
!      do k = 1,0 !nunitc
!      do j = 1,0 !nunitc
!      do i = 1,0 !nunitc
!      do m=1,natom_unit
!         natom = natom + 1
!         atom(natom) = atom(m)
!         rxyz(natom,1) = rxyz(m,1)+(i-1)*CL
!         rxyz(natom,2) = rxyz(m,2)+(j-1)*CL
!         rxyz(natom,3) = rxyz(m,3)+(k-1)*CL
!      end do
!      end do
!      end do
!      end do
!      print*, 'natom = ', natom
! origin at centre of box
!      rxyz(1:natom,1:3)=rxyz(1:natom,1:3) - CL*0.5_wp*(nunitc-1)
!      forall(i=1:natom) rxyz(i,1:3)=rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)
!      call write_xyz(14,0)
!      close(14,status='KEEP')

! determine connectivity
      do k=1,natom
         n=0
         do j=1,natom
            if(j == k) cycle
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.1_wp*bondl_SiO)**2) then
            !if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
                n=n+1
                proximity(k,n)=j
            end if
         end do
      end do
      do i=1,natom
         write (15,*) atom_name(atom(i)),proximity(i,1:4)
      end do
      close(15,status='KEEP')
      
      call new_frame_file(imve,'frame',0)
      call write_frame(imve,1,natom)



      deallocate(rxyz)
      deallocate(atom)
      allocate(rxyz(natomc*27,3),atom(natomc*27))
      rxyz(1:natomc,1:3) = rxyzc(1:natomc,1:3)
      atom(1:natomc) = atomc(1:natomc)


      natom = natomc
      proximity = 0
      do k = 1,natom
         n = 0
         do j=1,natom
            if(j == k) cycle
            if(atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.1_wp*bondl_SiO)**2) then 
           ! if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
               n = n+1
               proximity(k,n) = j
            end if
         end do
      end do
      call new_frame_file(imve,'frame',2)
      call write_frame(imve,1,natom)
!
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster



      natom = natomc
      proximity = 0
      do k = 1,natom
         n = 0
         do j=1,natom
            if(j == k) cycle
            if(atom(k) == atom(j)) cycle
            r3(:) = rxyz(j,:) - rxyz(k,:)
            !r3(:) = r3(:) - boxl*anint(r3(:)*boxli)
            CALL pbc(r3)
            ! MAY NEED TO CHANGE THE MULTIPLE OF BONDL_SiO
            if (dot_product(r3,r3) < (1.1_wp*bondl_SiO)**2) then
           ! if ( (len_3d(r3) <= bl1).OR.(len_3d(r3) <= bl2) ) then
               n = n+1
               proximity(k,n) = j
            end if
         end do
      end do
      call new_frame_file(imve,'frame',2)
      call write_frame(imve,1,natom)




      ib = 0
      print *,'natomc = ',natomc 
      do iz = -1,1
      do iy = -1,1
      do ix = -1,1
         if((ix==0).and.(iy==0).and.(iz==0))cycle
         ib = ib + 1
print *,ib,ib*natomc+1,(ib+1)*natomc
         rxyz(ib*natomc+1:(ib+1)*natomc,1) = rxyzc(1:natomc,1) + boxl*ix
         rxyz(ib*natomc+1:(ib+1)*natomc,2) = rxyzc(1:natomc,2) + boxl*iy
         rxyz(ib*natomc+1:(ib+1)*natomc,3) = rxyzc(1:natomc,3) + boxl*iz
         atom(ib*natomc+1:(ib+1)*natomc) = atomc(1:natomc)
      end do
      end do
      end do
      natom = 27*natomc
      print *,'natom = ',natom



!CONTAINS
!
!   pure function len_3d(r)
!      real(wp):: len_3d
!      real(wp),intent(in):: r(3)
!        len_3d = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
!   end function

END PROGRAM SIMBOX

