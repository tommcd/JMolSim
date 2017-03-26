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
include 'nlist_mod.f90'
!include 'crossbond_mod_v4.f90'
include 'frames_mod.f90'


PROGRAM crossbond
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
   USE nlist_mod
!   USE crossbond_mod
   USE files_mod
   IMPLICIT NONE
   type(atom_list):: lst_positive_x,lst_negative_x
   type(atom_list):: lst_positive_y,lst_negative_y
   type(atom_list):: lst_positive_z,lst_negative_z
   type(atom_list):: lst_pos_xy,lst_neg_xy
   type(atom_list):: lst_pos_xz,lst_neg_xz
   type(atom_list):: lst_pos_yz,lst_neg_yz
   type(atom_list):: lst_pos_xyz,lst_neg_xyz
   integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
   integer,allocatable:: special_xy1(:),special_xy2(:),special_xz1(:),special_xz2(:)
   integer,allocatable:: special_yz1(:),special_yz2(:),special_xyz1(:),special_xyz2(:)
   integer:: n_special_xy1,n_special_xy2,n_special_xz1,n_special_xz2
   integer:: n_special_yz1,n_special_yz2,n_special_xyz1,n_special_xyz2
   integer:: natomc
   integer,allocatable:: ibond_temp(:,:)
   integer,allocatable:: crossbond_x_temp(:),crossbond_y_temp(:),crossbond_z_temp(:)! new array definitions
   integer:: ineg_x,ineg_y,ineg_z,len_sq!,nbond_cell
   integer,allocatable:: lst_neg_x_face(:),lst_neg_y_face(:),lst_neg_z_face(:)
   integer,allocatable:: lst_pos_x_face(:),lst_pos_y_face(:),lst_pos_z_face(:)
   integer,allocatable:: lst_neg_spec_xy(:),lst_neg_spec_xz(:),lst_neg_spec_yz(:),lst_neg_spec_xyz(:)
   integer,allocatable:: lst_pos_spec_xy(:),lst_pos_spec_xz(:),lst_pos_spec_yz(:),lst_pos_spec_xyz(:)
   integer:: ipos_x,ipos_y,ipos_z,nbond_super_cell
   integer:: kk,n_crossbondtot,atomnum1,atomnum2,natom_copy,ifirst
   logical:: negx,negy,negz,ok,consistent! logical variables for neg faces
   integer,parameter:: natom_unit = 24
   integer,parameter:: nl = 1 
   integer,parameter:: nunitc(3) = 2*nl + 1
   integer:: n2ixyz((2*nl+1)**3,3),ni,ib_temp,ip,num_cell
   integer:: ixyz2n(-nl:nl,-nl:nl,-nl:nl) ! note:- nl should never be greater than 3
   real(wp):: timep(0:10)   
   real(wp):: r3(3),rj(3),bondl_orig=1.55_wp,bondl_scale
   real,allocatable:: rxyzc(:,:)
   integer,allocatable:: atomc(:)!,ia(:)
   real(wp):: CL,x,y,z,dr(3)!,bl1,bl2
   integer:: j,i,n,ib,k,imve,ix,iy,iz,iu3,ic,jj,nc,iu!,natom
   integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z
   character(32):: ctmp

   natom_max = natom_unit*(nunitc(1)*nunitc(2)*nunitc(3))**3
   allocate(rxyz(natom_max,3))
   allocate(atom(2*natom_max))
   
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
   
   allocate(proximity(natom_max,4))
!   allocate(atomc(1:natomc))
!   allocate(ia(10))
   proximity = 0
   nattached = 0
  
   CL = 7.17_wp/angstrom
   bondl_scale = bondl_SiO/bondl_orig
   CL = CL*bondl_scale
   boxl(:) = CL*nunitc(:)
   boxl2 = boxl/2.0_wp
   boxli = 1.0_wp/boxl
   
   num_cell = (2*nl+1)**3
   len_sq = (2*nl+1)**2

   allocate(lst_neg_x_face(len_sq),lst_neg_y_face(len_sq),lst_neg_z_face(len_sq))
   allocate(lst_pos_x_face(len_sq),lst_pos_y_face(len_sq),lst_pos_z_face(len_sq))
   ! 8 Si's in unit cell by number images by 4. because 4 bonds by Si 
   ! the number of bonds in the original imaged cell is as follows


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! will need to allocate prox. array to be a lot bigger than it is now
! rather than 3x3x3 it should be about 9x9x9

      CALL INIT_NLIST(boxl,5.0_wp/angstrom)

! image the original unit cell
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
      natom = natom_unit*nunitc(1)*nunitc(2)*nunitc(3)

      CALL NEW_NLIST
      
! the following is to get the connectivity using the nlist
      do i = 1, natom
         r3(1:3) = rxyz(i,1:3)
         ic = CELL(r3)
         CALL NEIGCELL(ic,2,neigh,ncell)
         n = 0
         cell_loop: do jj = 1, neigh
            nc = ncell(jj)
            if(nc==0) cycle cell_loop
            j=HOC(nc)
            cell_atom_loop: do while (j/=0)
               dr=r3-rxyz(j,1:3)
               if (j==i) GOTO 10
               ! if (dot_product(dr,dr)<=0.01)  cycle cell_loop
               CALL pbc(dr)
               if (dot_product(dr,dr)<(1.1*bondl_SiO)**2) then 
                  n=n+1
                  proximity(i,n) = j
                  if (n > 4) then
                     print*, 'n = ', n
                     print*, 'prox = ', proximity(i,:)
                  end if
               end if
10             CONTINUE
               j=LL(j)
            end do cell_atom_loop
         end do cell_loop
 
      end do 
      
      CALL bond_list

      call new_frame_file(iu,'frame',0)
      call write_frame(iu,1,natom)
      n_crossbondtot = 0.5_wp*nbondtot

      allocate(crossbond_x(nbondtot*num_cell),crossbond_y(nbondtot*num_cell))
      allocate(crossbond_z(nbondtot*num_cell))
      allocate(special_xy(nbondtot*num_cell),special_xz(nbondtot*num_cell))
      allocate(special_yz(nbondtot*num_cell),special_xyz(nbondtot*num_cell))
      CALL cpu_time(timep(0))
      ! Make arrays of crossing bonds
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i=1,nbondtot ! nbondtot calc'ed in bond_list
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
!need to have special cases for when the bonds cross the boundaries and in
!more than one plane. i.e. a crossbond in the x and y OR x and z OR y and z OR
!the case of crossing all three x and y and z
! special_xy,special_xz,special_yz,special_xyz
         
         if (ABS(r3(1)) > boxl2(1)) then
            if(ABS(r3(2)) > boxl2(2) .AND. ABS(r3(3)) > boxl2(3)) then
!!!! IGNORE FOR THE MOMENT AS THERE ARE NO INSTANCES OF THIS
               n_special_xyz = n_special_xyz + 1
               special_xyz(n_special_xyz) = i
               print*, 'special_xyz'
            else if (ABS(r3(2)) > boxl2(2)) then
               n_special_xy = n_special_xy + 1
               special_xy(n_special_xy) = i
            else if (ABS(r3(3)) > boxl2(3)) then
               n_special_xz = n_special_xz + 1
               special_xz(n_special_xz) = i
            else 
               n_crossbond_x = n_crossbond_x+1
               crossbond_x(n_crossbond_x)=i    
            end if
         else if (ABS(r3(2)) > boxl2(2)) then
            if (ABS(r3(3)) > boxl2(3)) then
               n_special_yz = n_special_yz + 1
               special_yz(n_special_yz) = i
            else 
               n_crossbond_y = n_crossbond_y+1
               crossbond_y(n_crossbond_y)=i  
            end if
         else if (ABS(r3(3)) > boxl2(3)) then
            n_crossbond_z = n_crossbond_z+1
            crossbond_z(n_crossbond_z)=i
         end if 
         
      end do
      
print*, 'special_xy = ',special_xy(1:20)
print*, 'special_xz = ',special_xz(1:20)
print*, 'special_yz = ',special_yz(1:20)
print*, 'special_xyz = ',special_xyz(1:20)
!stop
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
      do i = 1, n_special_xy
         ib = special_xy(i)
         if (rxyz(ibond(1,ib),1) > 0.0) then
            call add_to_list(lst_pos_xy,ibond(1,ib))
            call add_to_list(lst_neg_xy,ibond(2,ib))
            if(rxyz(ibond(2,ib),1) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               write(*,*) rxyz(ibond(2,ib),1),rxyz(ibond(1,ib),1)
               write(*,*) rxyz(ibond(2,ib),2),rxyz(ibond(1,ib),2)               
               stop 'error: bond does not cross boundary 1'
            end if
         else
            call add_to_list(lst_neg_xy,ibond(1,ib))
            call add_to_list(lst_pos_xy,ibond(2,ib))
         end if
      end do
!
      do i = 1, n_special_xz
         ib = special_xz(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_pos_xz,ibond(1,ib))
            call add_to_list(lst_neg_xz,ibond(2,ib))
            if(rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary 2'
            end if
         else
            call add_to_list(lst_neg_xz,ibond(1,ib))
            call add_to_list(lst_pos_xz,ibond(2,ib))
         end if
      end do
!
      do i = 1, n_special_yz
         ib = special_yz(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_pos_yz,ibond(1,ib))
            call add_to_list(lst_neg_yz,ibond(2,ib))
            if(rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary 3'
            end if
         else
            call add_to_list(lst_neg_yz,ibond(1,ib))
            call add_to_list(lst_pos_yz,ibond(2,ib))
         end if
      end do
!
      do i = 1, n_special_xyz
         ib = special_yz(i)
         if (rxyz(ibond(1,ib),3) > 0.0) then
            call add_to_list(lst_pos_xyz,ibond(1,ib))
            call add_to_list(lst_neg_xyz,ibond(2,ib))
            if(rxyz(ibond(2,ib),3) > 0.0) then
               write(*,*) rxyz(ibond(2,ib),3),rxyz(ibond(1,ib),3)
               stop 'error: bond does not cross boundary 4'
            end if
         else
            call add_to_list(lst_neg_xyz,ibond(1,ib))
            call add_to_list(lst_pos_xyz,ibond(2,ib))
         end if
      end do
         
!###############################################################
! the second imaging process takes place here imaging the 3x3x3 
! unit cell 
!###############################################################
! image the 3x3x3 box so it is essentially a 9x9x9 ~18,000 atoms
      ib = 0
      do iz = -nl,nl
      do iy = -nl,nl
      do ix = -nl,nl
         if((ix==0).and.(iy==0).and.(iz==0))cycle
         ib = ib + 1
         rxyz(ib*natom+1:(ib+1)*natom,1) = rxyz(1:natom,1) + boxl(1)*ix
         rxyz(ib*natom+1:(ib+1)*natom,2) = rxyz(1:natom,2) + boxl(2)*iy
         rxyz(ib*natom+1:(ib+1)*natom,3) = rxyz(1:natom,3) + boxl(3)*iz
         atom(ib*natom+1:(ib+1)*natom) = atom(1:natom)
      end do
      end do
      end do
      
      boxl(1) = boxl(1)*nunitc(1)
      boxl(2) = boxl(2)*nunitc(2)
      boxl(3) = boxl(3)*nunitc(3)
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
! the following line has been commented so that natom has the value of the !
! intermediate image (i.e. 3x3x3)                                          !
      natomc = natom*nunitc(1)*nunitc(2)*nunitc(3)                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
natom_copy = natom
natom = natomc
call get_free_file_unit(iu)
OPEN(UNIT=iu, FILE='imaged_cell.xyz')
call write_xyz(iu,0)
CLOSE(iu)
natom = natom_copy

! simply imaging the coordinates is insufficient. because of this the program
! will (should) image the proximity calculated for the 3x3x3 system as well as
! imaging the respective bond list and crossbond data. 
! As this imaging process continues conditional statements are used such that
! all the negative faces of the image are found and the cube numbers are stored.
! This will allow the appropriate crossbond correction to be applied to the
! bonds in each cube. 
! To aid this process two coupled arrays are filled during the imaging of data.
! These arrays are "ixyz2n", dim. (nunitc)x(nunitc)x(nunitc) and "n2ixyz", dim. 
! (nunitc**3)x3. Using these arrays it is possible to tell the atoms numbers
! used in a specific cell by entering the position relevant to original cell
! or by entering the cell number it can tell you the position relevant to
! original cell.
      ipos_x = 0;ipos_y =0;ipos_z =0;ineg_x = 0;ineg_y = 0;ineg_z = 0
      ib = 0
      ni = 1
      do iz = -nl,nl
      do iy = -nl,nl
      do ix = -nl,nl
         if((ix==0).and.(iy==0).and.(iz==0)) then
            n2ixyz(1,:) = (/0,0,0/)
            ixyz2n(ix,iy,iz) = 1
            cycle     
         end if
         ni = ni + 1
         ib = ib + 1         
         ixyz2n(ix,iy,iz) = ni
         n2ixyz(ni,:) = (/ix,iy,iz/)
!image the proximity array also
         do ip = 1, natom
            proximity(ib*natom+ip,1) = proximity(ip,1)+(ib*natom)
            proximity(ib*natom+ip,2) = proximity(ip,2)+(ib*natom)
            if (proximity(ip,3)/=0) then         
            proximity(ib*natom+ip,3) = proximity(ip,3)+(ib*natom)
            end if
            if (proximity(ip,4)/=0) then
            proximity(ib*natom+ip,4) = proximity(ip,4)+(ib*natom)
            end if
         end do  
         
if(ib==1) then
  call get_free_file_unit(iu)
  open(unit=iu,file='prox_data1.out')
  do i = 1, natomc
  write(iu,'(4(i6,2x))') proximity(i,:)
  end do
  close(iu)
end if

         ibond(:,ib*nbondtot+1:(ib+1)*nbondtot) = ibond(:,1:nbondtot)+ib*natom

! crossbond deals with actual bonds so a multiple of nbondtot must be added when
! imaging
         crossbond_x(ib*nbondtot+1:(ib+1)*nbondtot) = &
                                    crossbond_x(1:nbondtot) + ib*nbondtot
         crossbond_y(ib*nbondtot+1:(ib+1)*nbondtot) = &
                                    crossbond_y(1:nbondtot) + ib*nbondtot  
         crossbond_z(ib*nbondtot+1:(ib+1)*nbondtot) = &
                                    crossbond_z(1:nbondtot) + ib*nbondtot 
                                    
! the lst_neg/pos_x/y/z%ind are related to atom number so a multiple of natom
! must be added when imaging
         lst_negative_x%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_negative_x%ind(1:nbondtot) + ib*natom
         lst_negative_y%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_negative_y%ind(1:nbondtot) + ib*natom
         lst_negative_z%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_negative_z%ind(1:nbondtot) + ib*natom

         lst_positive_x%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_positive_x%ind(1:nbondtot) + ib*natom
         lst_positive_y%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_positive_y%ind(1:nbondtot) + ib*natom
         lst_positive_z%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_positive_z%ind(1:nbondtot) + ib*natom
                            
! need to image the special case atoms that cross more than one plane
         lst_neg_xy%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_neg_xy%ind(1:nbondtot) + ib*natom
         lst_neg_xz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_neg_xz%ind(1:nbondtot) + ib*natom
         lst_neg_yz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_neg_yz%ind(1:nbondtot) + ib*natom   
         lst_neg_xyz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_neg_xyz%ind(1:nbondtot) + ib*natom     

         lst_pos_xy%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_pos_xy%ind(1:nbondtot) + ib*natom
         lst_pos_xz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_pos_xz%ind(1:nbondtot) + ib*natom
         lst_pos_yz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_pos_yz%ind(1:nbondtot) + ib*natom
         lst_pos_xyz%ind(ib*nbondtot+1:(ib+1)*nbondtot) = &
                               lst_pos_xyz%ind(1:nbondtot) + ib*natom                               
                               
!          CONDITION A
         if(iy==-nl) then
            ineg_y = ineg_y + 1
            lst_neg_y_face(ineg_y) = ni
            
         else if(iy==-nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ni
            
         else if(iy==-nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ni
         endif
         
         
!          CONDITION B
         if(iy>=-nl .AND. iy<=nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ni
            
         else if(iy>=-nl .AND. iy<=nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ni 
         end if
         
         
!          CONDITION C
         if (iy==nl) then
            ipos_y = ipos_y + 1
            lst_pos_y_face(ipos_y) = ni
            
         else if (iy==nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ni
            
         else if (iy==nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ni
         end if
         
         if (iz==-nl) then
            ineg_z = ineg_z + 1
            lst_neg_z_face(ineg_z) = ni
         
         else if (iz==nl) then
            ipos_z = ipos_z + 1
            lst_pos_z_face(ipos_z) = ni
         end if
                                   
      end do
      end do
      end do


     
! the following loops go over all of the cells reassiging bonds that cross into
! neighbouring cells. the assumption is that all crossbonds can be defined in
! one of the x, y or z directions and no more. using this assumption the
! following uses the negative set of the bonds and reassigns the proximity
! arrays of the atoms. it treats the cubes with outer faces differently than
! those that have more than one internal edge. this means that there is a loop
! to deal with bonds going across the entire length of the cell and those that
! are bonded to physical neighbours.
print*, 'lst_neg_z_face = ',lst_neg_z_face
print*, 'natom = ', natom
print*, 'proximity(435,:) = ',  proximity(435,:)
                  
      do ib = 1, num_cell
         negx = ANY(lst_neg_x_face==ib)
         negy = ANY(lst_neg_y_face==ib)
         negz = ANY(lst_neg_z_face==ib)
!
!         print*, 'negz = ', negz
!print*, 'lst_negative_x%ind = ', lst_negative_x%ind((ib-1)*nbondtot)
!print*, 'lst_negative_x%ind = ', lst_negative_x%ind((ib-1)*nbondtot+1:(ib-1)*nbondtot+60) 
!print*, 'lst_negative_y%ind = ', lst_negative_y%ind((ib-1)*nbondtot+1:(ib-1)*nbondtot+60)
!print*, 'lst_negative_z%ind = ', lst_negative_z%ind((ib-1)*nbondtot+1:(ib-1)*nbondtot+60)

! treat the special crossers first
         if (negx) then
            if (negy) then
               do kk = 1, lst_neg_xy%n
                  atomnum1 = lst_neg_xy%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xy%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n((2*nl+n2ixyz(ib,1)),(2*nl+n2ixyz(ib,2)),n2ixyz(ib,3))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do
            else 
               do kk = 1, lst_neg_xy%n
                  atomnum1 = lst_neg_xy%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xy%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n((2*nl+n2ixyz(ib,1)),n2ixyz(ib,2)-1,n2ixyz(ib,3))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do               
            end if
            if (negz) then
               do kk = 1, lst_neg_xz%n
                  atomnum1 = lst_neg_xz%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xz%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n((2*nl+n2ixyz(ib,1)),n2ixyz(ib,2),(2*nl+n2ixyz(ib,3)))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do
            else
               do kk = 1, lst_neg_xz%n
                  atomnum1 = lst_neg_xz%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xz%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n((2*nl+n2ixyz(ib,1)),n2ixyz(ib,2),n2ixyz(ib,3)-1)
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do               
            end if               
         else
            if (negy) then
               do kk = 1, lst_neg_xy%n
                  atomnum1 = lst_neg_xy%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xy%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n(n2ixyz(ib,1)-1,(2*nl+n2ixyz(ib,2)),n2ixyz(ib,3))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do
            else 
               do kk = 1, lst_neg_xy%n
                  atomnum1 = lst_neg_xy%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xy%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n(n2ixyz(ib,1)-1,n2ixyz(ib,2)-1,n2ixyz(ib,3))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do                
            end if               
            if (negz) then
               do kk = 1, lst_neg_xz%n
                  atomnum1 = lst_neg_xz%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xz%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n(n2ixyz(ib,1)-1,n2ixyz(ib,2),(2*nl+n2ixyz(ib,3)))
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do
            else
               do kk = 1, lst_neg_xz%n
                  atomnum1 = lst_neg_xz%ind((ib-1)*nbondtot+kk)
                  atomnum2 = lst_pos_xz%ind((ib-1)*nbondtot+kk)
                  if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
                  nc = ixyz2n(n2ixyz(ib,1)-1,n2ixyz(ib,2),n2ixyz(ib,3)-1)
                  CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
                  CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)                     
               end do               
            end if  
         end if
      


!       x_direction 
!       first for the face cubes which involve the large crossbonds
         if(negx) then
            do kk = 1, lst_negative_x%n
               atomnum1 = lst_negative_x%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_x%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE               
              ! if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n((2*nl+n2ixyz(ib,1)),n2ixyz(ib,2),n2ixyz(ib,3))
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)    
               do i = 1,lst_negative_y%n
                  if (atomnum1==lst_negative_y%ind((ib-1)*nbondtot+i)  &
                     .AND.atomnum2==lst_positive_y%ind((ib-1)*nbondtot+i)) then
                     lst_negative_y%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_y%ind((ib-1)*nbondtot+i) = 0
                  elseif (atomnum2==lst_negative_y%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum1==lst_positive_y%ind((ib-1)*nbondtot+i)) then
                     lst_negative_y%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_y%ind((ib-1)*nbondtot+i) = 0
                  end if
               end do 
               do i = 1,lst_negative_z%n
                  if (atomnum1==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum2==lst_positive_z%ind((ib-1)*nbondtot+i)) then
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0
                  elseif (atomnum2==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum1==lst_positive_z%ind((ib-1)*nbondtot+i)) then     
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0
                  end if
               end do  
            end do
         elseif (.NOT. negx) then
!       Now for the inner boxes consistent that involve the crossbonds with a subsequent and 
!       preceeding box. 
            do kk = 1, lst_negative_x%n
               atomnum1 = lst_negative_x%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_x%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
               !if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n(n2ixyz(ib,1)-1,n2ixyz(ib,2),n2ixyz(ib,3))
!print*, 'b proximity(',atomnum1,':)',proximity(atomnum1,:) 
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom)
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)
!print*, 'a proximity(',atomnum1,':)',proximity(atomnum1,:) 
!print*, 'proximity(',atomnum2,':)',proximity(atomnum2,:)
!print*, 'proximity(',atomnum2+(nc-ib)*natom,':)',proximity(atomnum2+(nc-ib)*natom,:) 
!print*, ''
               do i = 1, lst_negative_y%n
                  if (atomnum1==lst_negative_y%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum2==lst_positive_y%ind((ib-1)*nbondtot+i)) then
                     lst_negative_y%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_y%ind((ib-1)*nbondtot+i) = 0
                  elseif (atomnum2==lst_negative_y%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum1==lst_positive_y%ind((ib-1)*nbondtot+i)) then
                     lst_negative_y%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_y%ind((ib-1)*nbondtot+i) = 0
                  end if
               end do
               do i = 1, lst_negative_z%n
                  if (atomnum1==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum2==lst_positive_z%ind((ib-1)*nbondtot+i)) then
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0                     
                  elseif (atomnum2==lst_negative_z%ind((ib-1)*nbondtot+i) & 
                     .AND.atomnum1==lst_positive_z%ind((ib-1)*nbondtot+i)) then   
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0 
                  end if
               end do   
            end do
         end if

if (ANY(proximity < 0)) stop 'stop 1' 

!       y_direction
         if(negy) then
            do kk = 1, lst_negative_y%n
               atomnum1 = lst_negative_y%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_y%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
               !if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n(n2ixyz(ib,1),2*nl+n2ixyz(ib,2),n2ixyz(ib,3))
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom)
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)
               do i = 1, lst_negative_z%n
                  if (atomnum1==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum2==lst_positive_z%ind((ib-1)*nbondtot+i)) then
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0 
                  elseif (atomnum2==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum1==lst_positive_z%ind((ib-1)*nbondtot+i)) then     
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0 
                  end if
               end do            
            end do
         elseif (.NOT. negy) then    
!       Now for the inner boxes that involve the crossbonds with a subsequent and 
!       preceeding box. 
            do kk = 1,  lst_negative_y%n
               atomnum1 = lst_negative_y%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_y%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
               !if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n(n2ixyz(ib,1),n2ixyz(ib,2)-1,n2ixyz(ib,3))
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom)
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)
               do i = 1, lst_negative_z%n
                  if (atomnum1==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum2==lst_positive_z%ind((ib-1)*nbondtot+i)) then
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0 
                  elseif (atomnum2==lst_negative_z%ind((ib-1)*nbondtot+i) &
                     .AND.atomnum1==lst_positive_z%ind((ib-1)*nbondtot+i)) then     
                     lst_negative_z%ind((ib-1)*nbondtot+i) = 0
                     lst_positive_z%ind((ib-1)*nbondtot+i) = 0 
                  end if
               end do               
            end do
         end if
         
if (ANY(proximity < 0)) stop 'stop 2' 

!       z_direction        
         if(negz) then  
!print*, 'lst_negative_z%ind(',ib,') = ', lst_negative_z%ind((ib-1)*nbondtot+1:(ib-1)*nbondtot+50)
            do kk = 1, lst_negative_z%n
               atomnum1 = lst_negative_z%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_z%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
!               if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n(n2ixyz(ib,1),n2ixyz(ib,2),2*nl+n2ixyz(ib,3))               
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom)
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)  
            end do
         else !if (negz) then      
!       Now for the inner boxes that involve the crossbonds with a subsequent and 
!       preceeding box.
!print*, 'ib = ', ib
!print*, 'lst_negative_z%ind(',ib,') = ', lst_negative_z%ind((ib-1)*nbondtot+1:(ib-1)*nbondtot+50)
            do kk = 1, lst_negative_z%n
               atomnum1 = lst_negative_z%ind((ib-1)*nbondtot+kk)
               atomnum2 = lst_positive_z%ind((ib-1)*nbondtot+kk)
               if (atomnum1 == 0 .OR. atomnum2 == 0) CYCLE
               !if(atomnum1==0 .AND. atomnum2==0) EXIT
               nc = ixyz2n(n2ixyz(ib,1),n2ixyz(ib,2),n2ixyz(ib,3)-1)
!if (ib == 23) then
!   print*, 'ib = ', ib
!   print*, 'atomnum1',atomnum1   
!   print*, 'atomnum2',atomnum2
!   print*, 'atomnum2+(',nc,'-',ib,')*natom = ',atomnum2+(nc-ib)*natom       
!   print*, 'proximity(',atomnum1,',:) = ',  proximity(atomnum1,:)      
!   print*, 'proximity(',atomnum2,',:) = ',  proximity(atomnum2,:) 
!   print*, 'proximity(',atomnum2+(nc-ib)*natom ,',:) = ',  proximity(atomnum2+(nc-ib)*natom ,:)      
!end if
!
!if (atomnum1 == 435) then
!   print*, 'ib = ', ib        
!   print*, 'atomnum1'
!   print*, 'proximity(435,:) = ',  proximity(435,:)
!   print*, ''     
!else if (atomnum2 == 435) then
!   print*, 'ib = ', ib        
!   print*, 'atomnum2',atomnum2
!   print*, 'atomnum1',atomnum1
!   print*, 'proximity(',atomnum2,',:) = ',  proximity(atomnum2,:)  
!   print*, 'proximity(',atomnum1,',:) = ',  proximity(atomnum1,:)    
!   print*, 'proximity(',atomnum2+(nc-ib)*natom ,',:) = ',  proximity(atomnum2+(nc-ib)*natom ,:)                     
!   print*, 'atomnum2+(',nc,'-',ib,')*natom = ',atomnum2+(nc-ib)*natom    
!   print*, ''            
!end if                
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(nc-ib)*natom) 
               CALL set_proximity(atomnum2+(nc-ib)*natom,atomnum1+(nc-ib)*natom,atomnum1)


            end do
         end if      
      end do    
!print*, 'proximity(435,:) = ',  proximity(435,:)  
!print*, 'atom 435 = ',atom_name(atom(435))
                  print*, ''  
      natom = natom*nunitc(1)*nunitc(2)*nunitc(3)      
      CALL cpu_time(timep(1))
      print*, 'time taken = ', timep(1) - timep(0)
      
      call check_proximity(consistent,ifirst,k)
      if (.NOT. consistent) then
         print*, 'consistent 5 = ',consistent
         print*, 'ifirst 5     = ',ifirst
         print*, 'proximity(',ifirst,':)',proximity(ifirst,:)
         print*, 'proximity(',k,':)',proximity(k,:)
         print*, 'proximity(',14353,':)',proximity(14353,:)
         !stop ': 5 proximity not consistent. fan...fucking...tastic!!!'
      end if
      call check_proximity2(consistent,ifirst)
      if (.NOT. consistent) then      
         print*, 'consistent 6 = ',consistent
         print*, 'ifirst 6     = ',ifirst   
         !stop ': 6 proximity not consistent. fan...fucking...tastic!!!'
      end if
      
      call bond_list
      call new_frame_file(iu,'done_and_dusted',1)
      call write_frame(iu,1,natom)
      
      
END PROGRAM crossbond 
