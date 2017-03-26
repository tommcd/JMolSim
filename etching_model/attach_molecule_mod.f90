MODULE attach_molecule_mod
   USE precision_mod
   USE global_vars_mod
   USE constants_mod
   USE atom_types_mod
   USE coordinates_mod
   USE connectivity_mod
   USE rand_mod
   implicit none

CONTAINS

   SUBROUTINE attach_mol(L,ndel,idel,bondl,bondang,natt,rmol,conect,atype,success)
!  SUBROUTINE attach_mol(L,bondl,bondang,pmol,success)
      USE nlist_mod
      USE rotate_axis_mod
      USE verlet_list_mod
      USE HKNonLattice_mod
      integer,intent(inout):: L
      integer,intent(in):: ndel,idel(:)
      real(wp),intent(in):: bondang,bondl
!     type(probe_mol):: pmol
      integer,intent(in):: natt
      real(wp),intent(in):: rmol(:,:)
      integer,intent(in):: conect(:,:)
      integer,intent(in):: atype(:)
      logical,intent(out):: success
!
      real(wp):: xx,yy,zz,phi,cphi,sphi,cosang,sinang,r2
      real(wp):: r3v(3),r12(3),r23(3),rSi(3),aa(3,3)
      integer:: i,j,LSi,ii,jj,M,ic,nc,nl,natom_old
      integer:: copproximity(natom_max,size(proximity,dim=2)),copatom(natom_max)
      real(wp):: coprxyz(natom_max,3)
      integer:: neighbors(343),nn
      logical:: overlap
!
      cosang = cos(bondang)
      sinang = sin(bondang)
      success = .false.
! save all the old data before the attempted attachment
      call NEW_VLIST
      do i = 1,natom
         coprxyz(i,1:3) = rxyz(i,1:3)
         copproximity(i,:) = proximity(i,:)
         copatom(i) = atom(i)
         cnlist(i) = nlist(i)
         clist(i,1:nlist(i)) = list(i,1:nlist(i))
      end do
      natom_old = natom
!
!-----Attaching a new molecule
      if (L > nseed) then
         do j = 1,ncmax(atom(L))
            ii = proximity(L,j)
            if ( atom(ii) == iSilicon ) then
               LSi = ii
               exit
            end if
         end do
         if ( j > ncmax(atom(L)) ) stop 'error: j > ncmax(atom(L)'
         r3v(1:3) = rxyz(L,1:3) - rxyz(LSi,1:3)
         call pbc(r3v)
         r12 = r3v/sqrt(dot_product(r3v,r3v))
      else ! attach to a seed
         r12 = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      end if

      do i = 1,ndel
         call delete_atom(idel(i))
      end do
      if ( L > natom) L = natom

!     prepare the position for silica
      call zaxis2vect(r12, aa)
      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      zz = -bondl*cosang
      xx = bondl*sinang*cphi
      yy = bondl*sinang*sphi
      r23 = matmul(aa,(/xx,yy,zz/))
      rSi(1:3) = rxyz(L,1:3) + r23
      r23 = r23/sqrt(dot_product(r23,r23))

!     prepare the position of the molecule
      call zaxis2vect(r23,aa)
!
!     attach the molecule(mol)
      phi = rand()*2.0_wp*pi
      cphi = cos(phi)
      sphi = sin(phi)
      do i = 1,natt
         rxyz(natom + i,1) = cphi*rmol(i,1) - sphi*rmol(i,2)
         rxyz(natom + i,2) = sphi*rmol(i,1) + cphi*rmol(i,2)
         rxyz(natom + i,3) = rmol(i,3)
      end do
      do i = 1,natt
         rxyz(natom + i,1:3) = matmul(aa,rxyz(natom + i,1:3))
      end do
      do i = 1,natt
         rxyz(natom + i,1) = rxyz(natom + i,1) + rSi(1)
         rxyz(natom + i,2) = rxyz(natom + i,2) + rSi(2)
         rxyz(natom + i,3) = rxyz(natom + i,3) + rSi(3)
      end do
!
!-----Apply PBC to the new atoms
      call pbc(rxyz(natom + 1:natom + natt,:))
!
      proximity(natom + 1:natom + natt,:) = conect(1:natt,:)
      where (proximity(natom + 1:natom + natt,:) /= 0) &
             proximity(natom + 1:natom + natt,:) = &
             proximity(natom + 1:natom + natt,:) + natom
      where (proximity(natom + 1,:) > (natom + natt)) proximity(natom + 1,:) = 0
      atom(natom + 1:natom + natt) = atype(1:natt)

      call set_proximity(natom + 1,0,L)
      call set_proximity(L,0,natom + 1)

      atom(L) = iOxygen
      atom(natom + 1:natom + natt) = atype(1:natt)
!
      atomL(natom + 1:natom + natt) = atomL(L)  ! cluster label for new atoms
      natom = natom + natt
      call NEW_NLIST(1,natom)
!
! Check for overlap
      overlap = .false.
      outer_atom_loop: do ii = natom - natt + 1,natom
      nl = nlayers_ll(sigLJ(atom(ii)))
      if (nl == 0) cycle
      ic = CELL(rxyz(ii,:))  !link list
      call NEIGCELL(ic,nl,nn,neighbors)
      cell_loop: do jj = 1,nn
         nc = neighbors(jj)
         if (nc == 0) cycle cell_loop
         M = HOC(nc)
         do while (M /= 0)
            if ( M > natom - natt) GOTO 100
            if (atom(M) == iSilicon) GOTO 100
            if (M == ii) GOTO 100
            if (nearest_neighbor2(ii,m)) GOTO 100
            r3v(1:3) = rxyz(M,1:3) - rxyz(ii,1:3)
            call pbc(r3v)
            r2 = r3v(1)**2 + r3v(2)**2 + r3v(3)**2
            if (r2 <= (sigLJ_2(atom(m)) + sigLJ_2(atom(ii)))**2) then
               overlap = .true.
               exit outer_atom_loop
            end if
100         M = LL(M)
         end do
      end do cell_loop
      end do outer_atom_loop
!
      if (.not.overlap) then
         success = .true.
         RETURN
      end if
!
999   CONTINUE
      atom(natom + 1:natom + natt) = 0
      proximity(natom + 1:natom + natt,:) = 0
      natom = natom_old
      do i = 1,natom
         rxyz(i,1:3) = coprxyz(i,1:3)
         proximity(i,:) = copproximity(i,:)
         atom(i) = copatom(i)
         nlist(i) = cnlist(i)
         list(i,1:cnlist(i)) = clist(i,1:cnlist(i))
      end do
      call NEW_NLIST(1,natom)
      RETURN
   END SUBROUTINE attach_mol

END MODULE attach_molecule_mod

