
MODULE insert_mod
    USE precision_mod, only: wp
    USE global_vars
    USE nlist_mod
    USE atom_types
    USE coordinates_mod
    USE seaton_mod
    USE rand_mod
    USE quat2mat_mod
    USE ran_point_sphere

    implicit none

CONTAINS

   SUBROUTINE insert_mol_N2_O2(n_mol,real_xyz,ntrial,rad1)

      integer,intent(in):: ntrial
      real(wp) fover
      real(wp):: aa(3,3),q(4),x,y,z,real_xyz(:,:,:),rad1(3)
      real(wp):: mol_vec1(3,3)
      integer:: it,nat,i,n_mol,n
      

      mol_loop: do n = 1, n_mol
      
      insert_loop: do it = 1,ntrial
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         z = ( 1.0_wp - 2.0_wp*rand() )*boxl2

            !call random_rotation(aa)
            call ran4sph(q)
            call quat2mat(q,aa)
            do i = 1,3
            mol_vec1(:,i) = matmul(aa,mol_vec(:,i))
            end do
            
         mol_vec1(1,:) = mol_vec1(1,:) + x
         mol_vec1(2,:) = mol_vec1(2,:) + y
         mol_vec1(3,:) = mol_vec1(3,:) + z
!!       pr1%r(1:2,:) = pr1%r(1:2,:) - anint(pr1%r(1:2,:)*boxli)*boxl
         do i = 1,3
            call pbc(mol_vec1(1:3,i))
         end do
! first check if there is a large overlap
         
         if(.not.overlap(1.12_wp,rad1,mol_vec1)) then
          real_xyz(n,1,1) = mol_vec1(1,1)
          real_xyz(n,1,2) = mol_vec1(2,1)
          real_xyz(n,1,3) = mol_vec1(3,1)
          real_xyz(n,2,1) = mol_vec1(1,2)
          real_xyz(n,2,2) = mol_vec1(2,2)
          real_xyz(n,2,3) = mol_vec1(3,2)
     write(*,*) n     
                 cycle mol_loop
         end if
         
      end do insert_loop      
 end do mol_loop
      RETURN
   END SUBROUTINE insert_mol_N2_O2
   
   SUBROUTINE insert_mol_CO2(n_mol,real_xyz,ntrial,rad1)
      integer,intent(in):: ntrial
      real(wp):: fover
      real(wp):: aa(3,3),q(4),x,y,z,real_xyz(:,:,:),rad1(3)
      real(wp):: mol_vec1(3,3)
      integer:: it,nat,i,n_mol,n

      mol_loop: do n = 1, n_mol
      insert_loop: do it = 1,ntrial
         x = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         y = ( 1.0_wp - 2.0_wp*rand() )*boxl2
         z = ( 1.0_wp - 2.0_wp*rand() )*boxl2

            !call random_rotation(aa)
            call ran4sph(q)
            call quat2mat(q,aa)
            do i = 1,3
            mol_vec1(:,i) = matmul(aa,mol_vec(:,i))
            end do
            
         mol_vec1(1,:) = mol_vec1(1,:) + x
         mol_vec1(2,:) = mol_vec1(2,:) + y
         mol_vec1(3,:) = mol_vec1(3,:) + z
!!       pr1%r(1:2,:) = pr1%r(1:2,:) - anint(pr1%r(1:2,:)*boxli)*boxl
         do i = 1,3
            call pbc(mol_vec1(1:3,i))
         end do
! first check if there is a large overlap
         
         if(.not.overlap(1.12_wp,rad1,mol_vec1)) then
          real_xyz(n,1,1) = mol_vec1(1,1)
          real_xyz(n,1,2) = mol_vec1(2,1)
          real_xyz(n,1,3) = mol_vec1(3,1)
          real_xyz(n,2,1) = mol_vec1(1,2)
          real_xyz(n,2,2) = mol_vec1(2,2)
          real_xyz(n,2,3) = mol_vec1(3,2)         
          real_xyz(n,3,1) = mol_vec1(1,3)
          real_xyz(n,3,2) = mol_vec1(2,3)
          real_xyz(n,3,3) = mol_vec1(3,3) 
   write(*,*) n       
                 cycle mol_loop
         end if
         
      end do insert_loop      
 end do mol_loop                                           
      RETURN
!
   END SUBROUTINE insert_mol_CO2
!
       function overlap(fr,rad,mol_vec1)
         logical:: overlap
         integer:: k,jj,j,ic,nc 
         real(wp):: dr(3),rad(3),fr,mol_vec1(3,3)
         overlap = .false.
         probe_atom_loop: do k = 1,3
         ic = CELL( mol_vec1(1:3,k) )
         call NEIGCELL(ic,1,neigh,ncell)
         cell_loop: do jj = 1,neigh
            nc = ncell(jj)
            if (nc == 0) cycle cell_loop
            j = HOC(nc)
            cell_atom_loop: do while (j /= 0)
               if (atom(j) /= iSilicon) then
               dr(1:3) = mol_vec1(1:3,k) - rxyz(j,1:3)
               call pbc(dr)
               if (dot_product(dr,dr) < ( fr*(rad(k) + sigi2(atom(j)) ) )**2) then
                  overlap = .true.
                  return
               end if
               end if
               j = LL(j)
            end do cell_atom_loop
         end do cell_loop
         end do probe_atom_loop
      end function

       function overlap2(fr,rad,mol_vec1)
         logical:: overlap2
         integer:: k,jj,j,ic,nc 
         real(wp):: dr(3),rad(3),fr,mol_vec1(3,3)
         overlap2 = .false.
         
         probe_atom_loop: do k = 1,3
         
         do j = 1,natom
            if (atom(j) == iSilicon) cycle
            dr(1:3) = mol_vec1(1:3,k) - rxyz(j,1:3)
            call pbc(dr)
            if (dot_product(dr,dr) < ( fr*(rad(k) + sigi2(atom(j)) ) )**2) then
               overlap2 = .true.
               return
            end if
         end do
         
         end do probe_atom_loop
      end function
END MODULE insert_mod



