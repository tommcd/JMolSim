MODULE crossbond_mod
   USE precision_mod
   IMPLICIT NONE


CONTAINS   

   SUBROUTINE crossbond_correct(nl,crossbond_x,crossbond_y,crossbond_z)
      USE atom_list_mod
      USE bond_list_mod
      USE global_vars
      USE atom_types
      USE seaton_mod
      USE constants_mod
      USE connectivity_mod
      integer,allocatable:: ibond_temp(:,:)
      integer,allocatable:: crossbond_x_temp(:),crossbond_y_temp(:),crossbond_z_temp(:)! new array definitions
      integer:: ineg_x,ineg_y,ineg_z,len_sq,nbond_cell,natomc
      integer,allocatable:: lst_neg_x_face(:),lst_neg_y_face(:),lst_neg_z_face(:)
      integer,allocatable:: lst_pos_x_face(:),lst_pos_y_face(:),lst_pos_z_face(:)
      integer,allocatable:: lst_neg_x_temp(:),lst_neg_y_temp(:),lst_neg_z_temp(:)
      integer,allocatable:: lst_pos_x_temp(:),lst_pos_y_temp(:),lst_pos_z_temp(:)      
      integer,allocatable:: n_faces(:)
      integer:: ncell, ipos_x,ipos_y,ipos_z,nl,nbond_super_cell,ib,ix,iy,iz
      integer:: kk,n_crossbondtot,atomnum1,atomnum2
      logical:: negx,negy,negz ! logical variables for neg faces
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)

      type(atom_list):: lst_positive_x,lst_negative_x
      type(atom_list):: lst_positive_y,lst_negative_y
      type(atom_list):: lst_positive_z,lst_negative_z
      ncell = (2*ncell+1)**3
      len_sq = (2*nl+1)**2
      allocate(n_faces((2*nl+1)**3))
      allocate(lst_neg_x_face(len_sq),lst_neg_y_face(len_sq),lst_neg_z_face(len_sq))
      allocate(lst_pos_x_face(len_sq),lst_pos_y_face(len_sq),lst_pos_z_face(len_sq))
      ! 8 Si's in unit cell by number images by 4. because 4 bonds by Si 
      ! the number of bonds in the original imaged cell is as follows
      nbond_cell = (8*(2*nl+1)**3)*4  
      n_crossbondtot = 1000 ! 1000 should cover the num of crossbonds in any dir.
      nbond_super_cell = n_crossbondtot*(2*nl+1)**3  

print*, 'n_crossbondtot =', n_crossbondtot
      
! tidying up and reassigning sizes to the arrays in use 
      allocate(ibond_temp(2,natom_max))
      allocate(crossbond_x_temp(natom_max),crossbond_y_temp(natom_max))
      allocate(crossbond_z_temp(natom_max))    
      allocate(lst_neg_x_temp(n_crossbondtot),lst_neg_y_temp(n_crossbondtot))
      allocate(lst_neg_z_temp(n_crossbondtot))
      allocate(lst_pos_x_temp(n_crossbondtot),lst_pos_y_temp(n_crossbondtot))
      allocate(lst_pos_z_temp(n_crossbondtot))      

print*, 'size(crossbond_x) =', size(crossbond_x)
      ibond_temp = ibond 
      crossbond_x_temp(:) = crossbond_x(:)
      crossbond_y_temp(:) = crossbond_y(:)
      crossbond_z_temp(:) = crossbond_z(:)
print*, 'lst_negative_x%i(3) =',lst_negative_x
      lst_neg_x_temp(:) = lst_negative_x%i(:)
      lst_neg_y_temp(:) = lst_negative_y%i(:)
      lst_neg_z_temp(:) = lst_negative_z%i(:)
      lst_pos_x_temp(:) = lst_positive_x%i(:)
      lst_pos_y_temp(:) = lst_positive_y%i(:)
      lst_pos_z_temp(:) = lst_positive_z%i(:)      
      deallocate(ibond,crossbond_x,crossbond_y,crossbond_z)
      deallocate(lst_negative_x%i,lst_negative_y%i,lst_negative_z%i)
      deallocate(lst_positive_x%i,lst_positive_y%i,lst_positive_z%i)      
      
      allocate( ibond(2,nbond_cell*((2*nl+1)**3)) )
      allocate(crossbond_x(nbond_super_cell),crossbond_y(nbond_super_cell))
      allocate(crossbond_z(nbond_super_cell))
      
      allocate(lst_negative_x%i(nbond_super_cell))
      allocate(lst_negative_y%i(nbond_super_cell))
      allocate(lst_negative_z%i(nbond_super_cell))   
      allocate(lst_positive_x%i(nbond_super_cell))
      allocate(lst_positive_y%i(nbond_super_cell))
      allocate(lst_positive_z%i(nbond_super_cell))
      ibond(:,1:nbond_cell*4) = ibond_temp(:,1:nbond_cell*4) 
      crossbond_x(1:natom_max) = crossbond_x_temp(1:natom_max)
      crossbond_y(1:natom_max) = crossbond_y_temp(1:natom_max)
      crossbond_z(1:natom_max) = crossbond_z_temp(1:natom_max)
      lst_negative_x%i(1:natom_max) = lst_neg_x_temp(1:natom_max)
      lst_negative_y%i(1:natom_max) = lst_neg_x_temp(1:natom_max)
      lst_negative_z%i(1:natom_max) = lst_neg_x_temp(1:natom_max)
      lst_positive_x%i(1:natom_max) = lst_pos_x_temp(1:natom_max)
      lst_positive_y%i(1:natom_max) = lst_pos_x_temp(1:natom_max)
      lst_positive_z%i(1:natom_max) = lst_pos_x_temp(1:natom_max)      
! initialis the variables for neg face list assignment     
      ineg_x = 0     
      ineg_y = 0
      ineg_z = 0
      ipos_x = 0
      ipos_y = 0
      ipos_z = 0
!!!!!!!!  ASSUMING THAT NATOMC = NUMBER ATOMS IN IMAGED CELL NOT UNIT CELL. 
!!!!!!!! MAY NEED 
      natomc = natom 
      
! as well as imaging the coordinates of original box this will image the list to
! the images as well as the crossbonds and the negative and positive bond lists
! also it will fill the arrays for the negative/positive faces (pos not used)
      ib = 0
      nl = 1
      do iz = -nl,nl
      do iy = -nl,nl
      do ix = -nl,nl
         if((ix==0).and.(iy==0).and.(iz==0))cycle
         ib = ib + 1
         atom(ib*natomc+1:(ib+1)*natomc) = atom(1:natomc)
! so need to update the bond list as well as the crossbond list as the imaging
! is being performed. 
         ibond(:,ib*nbond_cell+1:(ib+1)*nbond_cell) = ibond(:,1:nbond_cell)+natomc
         crossbond_x(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                                    crossbond_x(1:n_crossbondtot) + n_crossbondtot
         crossbond_y(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                                    crossbond_y(1:n_crossbondtot) + n_crossbondtot  
         crossbond_z(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                                    crossbond_z(1:n_crossbondtot) + n_crossbondtot 

         lst_negative_x%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_negative_x%i(1:n_crossbondtot) + n_crossbondtot
         lst_negative_y%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_negative_y%i(1:n_crossbondtot) + n_crossbondtot
         lst_negative_z%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_negative_y%i(1:n_crossbondtot) + n_crossbondtot

         lst_positive_x%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_positive_x%i(1:n_crossbondtot) + n_crossbondtot
         lst_positive_y%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_positive_y%i(1:n_crossbondtot) + n_crossbondtot
         lst_positive_z%i(ib*n_crossbondtot+1:(ib+1)*n_crossbondtot) = &
                               lst_positive_y%i(1:n_crossbondtot) + n_crossbondtot                              
                                    
!          CONDITION A
         if(iy==-nl) then
            ineg_y = ineg_y + 1
            lst_neg_y_face(ineg_y) = ib
            
         else if(iy==-nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ib
            
         else if(iy==-nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ib
         endif
         
         
!          CONDITION B
         if(iy>=-nl .AND. iy<=nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ib
            
         else if(iy>=-nl .AND. iy<=nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ib 
         end if
         
         
!          CONDITION C
         if (iy==nl) then
            ipos_y = ipos_y + 1
            lst_pos_y_face(ipos_y) = ib
            
         else if (iy==nl .AND. ix==-nl) then
            ineg_x = ineg_x + 1
            lst_neg_x_face(ineg_x) = ib
            
         else if (iy==nl .AND. ix==nl) then
            ipos_x = ipos_x + 1
            lst_pos_x_face(ipos_x) = ib
         end if
         
         if (iz==-nl) then
            ineg_z = ineg_z + 1
            lst_neg_z_face(ineg_z) = ib
         
         else if (iz==nl) then
            ipos_z = ipos_z + 1
            lst_pos_z_face(ipos_z) = ib
         end if
                                    
      end do
      end do
      end do

      

! the crossbonds need only be calculated for the first box. using the
! crossbond list then the elements in it could be increased to give the
! correct crossbond list for a particular box by multiplying the atom
! numbers by the box index number. 
      
! the following loops go over all of the cells reassiging bonds that cross into
! neighbouring cells. the assumption is that all crossbonds can be defined in
! one of the x, y or z directions and no more. using this assumption the
! following uses the negative set of the bonds and reassigns the proximity
! arrays of the atoms. it treats the cubes with outer faces differently than
! those that have more than one internal edge. this means that there is a loop
! to deal with bonds going across the entire length of the cell and those that
! are bonded to physical neighbours.
      
      
      
      
      do ib = 1, ncell
         negx = ANY(lst_neg_x_face==ib)
         negy = ANY(lst_neg_x_face==ib)
         negz = ANY(lst_neg_x_face==ib)
         

      
!       x_direction
!       first for the face cubes which involve the large crossbonds
         if(negx) then
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_x%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_x%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2+(2*nl*natomc)) 
            end do
         else
!       Now for the inner boxes that involve the crossbonds with a subsequent and 
!       preceeding box. 
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_x%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_x%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2-natomc) 
            end do
         end if
         
!       y_direction
         if(negy) then
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_y%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_y%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2+ &
                                        ((2*nl+1)**2*natomc - (2*nl+1)*natomc))
            end do
         else      
!       Now for the inner boxes that involve the crossbonds with a subsequent and 
!       preceeding box. 
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_y%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_y%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2-((2*nl+1)*natomc)) 
            end do
         end if
      
      
         if(negz) then  
!       z_direction
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_z%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_z%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2+ &
                                                          ((2*nl+1)**3)*natomc)
            end do
         else      
!       Now for the inner boxes that involve the crossbonds with a subsequent and 
!       preceeding box. 
            do kk = 1, n_crossbondtot
               atomnum1 = lst_negative_z%i(ib*n_crossbondtot+kk)
               atomnum2 = lst_positive_z%i(ib*n_crossbondtot+kk)
               if(atomnum1==0) EXIT
               CALL set_proximity(atomnum1,atomnum2,atomnum2-&
                                                          ((2*nl+1)**2)*natomc) 
            end do
         end if
      
      end do    

   END SUBROUTINE crossbond_correct
END MODULE crossbond_mod
