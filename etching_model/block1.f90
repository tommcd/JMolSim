MODULE BLOCK1
      USE precision_mod, only: wp
      implicit none
      integer,parameter:: nmx = 8000
      !real(wp),parameter:: del_rxn = 0.14_wp
      !real(wp),parameter:: e_activ = 0.75_wp
      real(wp),parameter:: del_rxn = 0.2_wp
      real(wp),parameter:: e_activ = 1.0_wp
      integer,parameter:: ioxygen = 0
      integer,parameter:: isilicon = 1
      integer,parameter:: ihydrogen = 2
      real(wp),parameter:: angstrom = 10.0_wp  ! Angstroms/nm
      real(wp),parameter:: erg_ev = 6.241457E+11_wp
      real(wp),parameter:: K_ev = 8.6173423E-5_wp
      real(wp),parameter:: qstar = 1.19999_wp
      character(2),parameter:: atom_name(0:1) = (/ 'O ','Si' /)
      real(wp),parameter:: sigma_Obridge = 2.7_wp/angstrom
      real(wp),parameter:: sigma_H2O = 3.166_wp/angstrom
      real(wp),parameter:: sigma_OH = 3.00_wp/angstrom
      real(wp),parameter:: sigma_O = 2.94_wp/angstrom
      real(wp),parameter:: sigma_Si = 4.20_wp/angstrom
      real(wp),parameter:: sigma_C = 2.70_wp/angstrom
      real(wp),parameter:: sigma_N = 3.296_wp/angstrom
      real(wp),parameter:: bondl_O2 = 1.169_wp/angstrom
      real(wp),parameter:: bondl_N2 = 1.097_wp/angstrom
      real(wp),parameter:: bondl_SiO4 = 1.62_wp/angstrom
      real(wp),parameter:: bondl_CO2 = 0.958_wp/angstrom
      real(wp),parameter:: r_tet = bondl_SiO4 + sigma_OH*0.5_wp
      real(wp),parameter:: sigma(0:1) = (/ sigma_OH,sigma_Si /)
      real(wp),parameter:: sigma_2(0:1) = sigma*0.5_wp
      integer,parameter:: ncmax(0:1) = (/2,4/)
      real(wp):: rxyz(nmx,3),boxl,boxli,boxl2
      integer:: proximity(nmx,4),atom(nmx),natom,nseed,nattached
      integer:: danglingbond(nmx),nrelax
!
   CONTAINS

   SUBROUTINE set_proximity(iat,iold,inew)
      integer,intent(in):: iat,iold,inew
      integer:: i
      do i = 1,ncmax(atom(iat))
         if (proximity(iat,i) == iold) then
            proximity(iat,i) = inew
            exit
         end if
         if (i == ncmax(atom(iat))) STOP 'set_proximity: ERROR'
      end do
   END SUBROUTINE set_proximity

   SUBROUTINE check_proximity(consistent,ifirst)
      logical,intent(out):: consistent
      integer,intent(out):: ifirst
      integer:: i,j,k
      consistent = .true.
      ifirst = 0
      do i = 1,natom
         do j = 1,ncmax(atom(i))
            k = proximity(i,j)
            if (k == 0) cycle
            if ( count(proximity(k,:) == i) /= 1) then
               consistent = .false.
               ifirst = i
               RETURN
            end if
         end do
      end do
   END SUBROUTINE

   SUBROUTINE write_proximity(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,j,k
      write(iu,*)'# nattached = ',nattached
      do i = 1,natom
         write(iu,'(i9,1x,i4)') atom(i),danglingbond(i)
      end do
      write(iu,'(/)')
      do i = 1,natom
         k = count(proximity(i,:) > 0)
         write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))')&
            i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE write_proximity

   SUBROUTINE read_proximity(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,j,it
      character(32):: ctmp
      read(iu,*) ctmp,ctmp,ctmp,it
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) atom(i),danglingbond(i)
      end do
      do i = 1,natom
         read(iu,*) it,ctmp,it,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      rewind(iu)
   END SUBROUTINE read_proximity

   SUBROUTINE delete_atom(i)
      integer,intent(in):: i
      integer:: j,k
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         where( proximity(k,:) == i ) proximity(k,:) = 0
      end do
      rxyz(i,1:3) = rxyz(natom,1:3)
      proximity(i,:) = proximity(natom,:)
      atom(i) = atom(natom)
      danglingbond(i) = danglingbond(natom)
      do j = 1,ncmax(atom(i))
         k = proximity(i,j)
         if (k == 0) cycle
         where( proximity(k,:) == natom ) proximity(k,:) = i
      end do
      rxyz(natom,1:3) = 0.0_wp
      proximity(natom,:) = 0
      danglingbond(natom) = 0
      atom(natom) = 0
      natom = natom - 1
   END SUBROUTINE delete_atom

   PURE FUNCTION nearest_neighbor2(L,M)
      logical:: nearest_neighbor2
      integer,intent(in):: L,M
      integer:: i,j,LP
      nearest_neighbor2 = .false.
      do i = 1,ncmax(atom(L))
         LP = proximity(L,i)
         if (LP == 0) cycle
         if (LP == M) then
            nearest_neighbor2 = .true.
            RETURN
         end if
         do j = 1,ncmax(atom(LP))
            if (proximity(LP,j) == M) then
               nearest_neighbor2 = .true.
               RETURN
            end if
         end do
      end do
   END FUNCTION nearest_neighbor2

   PURE FUNCTION is_dangling_group(iSi,iat)
      integer,intent(in):: iSi,iat
      logical:: is_dangling_group
      integer:: k,j
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         if (j == iat) cycle
         if (j == 0) cycle
         if (danglingbond(j) /= 1) then
            is_dangling_group = .FALSE.
            RETURN
         end if
      end do
      is_dangling_group = .TRUE.
   END FUNCTION

   SUBROUTINE delete_group(iSi)
!!   SUBROUTINE delete_group(iSi,iat)
      integer,intent(in):: iSi !!,iat
      integer:: k,j
      do k = 1,ncmax(atom(iSi))
         j = proximity(iSi,k)
         !!if (j == iat) cycle
         if (j == 0) cycle
         if (danglingbond(j) /= 1) then
            write(*,*)'delete_group: error'
            stop 'not a dangling group'
         end if
         call delete_atom(j)
      end do
      call delete_atom(iSi)
   END SUBROUTINE

!   SUBROUTINE undelete_atom(i) ! not finished
!      integer,intent(in):: i
!      integer:: j,k
!      natom = natom + 1
!      atom(natom) = atom(i)
!      danglingbond(natom) = danglingbond(i)
!      proximity(natom,:) = proximity(i,:)
!      rxyz(natom,1:3) = rxyz(i,1:3)
!      do j = 1,ncmax(atom(i))
!         k = proximity(i,j)
!         if (k == 0) cycle
!         where( proximity(k,:) == i ) proximity(k,:) =natom
!      end do
!      danglingbond(i) = danglingbond_old(i)
!      atom(i) = atom_old(i)
!      proximity(i,:) = proximity_old(i,:)
!      rxyz(i,1:3) = rxyz_old(i,1:3)
!      do j = 1,ncmax(atom(i))
!         k = proximity(i,j)
!         if (k == 0) cycle
!         proximity(k,:) = copproximity_old(k,:)
!      end do
!   END SUBROUTINE undelete_atom

   SUBROUTINE WRITE_XYZ(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i
      write(iu,*)'nattached = ',nattached
      write(iu,*) natom
      write(iu,*)'Silica'
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8),1x,i3)') atom_name(atom(i)),rxyz(i,1)*Angstrom &
         ,rxyz(i,2)*Angstrom,rxyz(i,3)*Angstrom,danglingbond(i)
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,it
      real(wp):: x,y,z
      character(32):: ctmp
      READ(iu,*) ctmp,ctmp,it
      if (it /= nattached) stop 'it /= nattached'
      READ(iu,*) natom
      READ(iu,*) ctmp
      do i = 1,natom
         READ(iu,*) ctmp,x,y,z,danglingbond(i)
         rxyz(i,1) = x/Angstrom
         rxyz(i,2) = y/Angstrom
         rxyz(i,3) = z/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE

