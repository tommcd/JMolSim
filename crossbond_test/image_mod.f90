MODULE image_mod
USE precision_mod
USE coordinates_mod
USE global_vars
USE frames_mod
IMPLICIT NONE

CONTAINS 
   
   SUBROUTINE image_box(boxl,rxyz)
   real(wp) :: boxl,boxl2
   real(wp) :: natom, rxyz(:,:)
   integer :: atom(:)
   boxl2 = boxl/2

      ! initial translation of box occurs as follows 
      rxyz(1:natom,1:2) = rxyz(1:natom,1:2) + boxl2

      ! image in y direction
      rxyz(natom+1:2*natom,2) = rxyz(1:natom,2) - boxl
      atom(natom+1:2*natom) = atom(1:natom)      
      
      ! image in x direction 
      rxyz(2*natom+1:4*natom,1) = rxyz(1:2*natom,1) - boxl
      atom(2*natom+1:4*natom) = atom(1:2*natom)
      
      ! image in z direction
      rxyz(4*natom+1:8*natom,3) = rxyz(1:4*natom,3) + boxl
      rxyz(8*natom+1:12*natom,3) = rxyz(1:4*natom,3) - boxl
      atom(4*natom+1:8*natom) = atom(1:4*natom)     
      atom(8*natom+1:12*natom) = atom(1:4*natom)
   

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



      
   END SUBROUTINE image_box
   
END MODULE image_mod

