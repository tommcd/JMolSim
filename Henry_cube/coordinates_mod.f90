

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:),Ox_xyz(:,:,:)
   real(wp),allocatable:: Ox_fxyz(:,:,:),N2_fxyz(:,:,:),CO2_fxyz(:,:,:)
   real(wp),allocatable:: Imige_Ox_xyz(:,:,:),N2_xyz(:,:,:),CO2_xyz(:,:,:)
   real(wp),allocatable:: Imige_CO2_xyz(:,:,:),Imige_N2_xyz(:,:,:) 
   real(wp),allocatable:: Ox_vxyz(:,:,:),N2_vxyz(:,:,:),CO2_vxyz(:,:,:),mol_vec(:,:)
   
   real(wp):: boxl,boxli,boxl2
   interface pbc
      module procedure  pbc_v, pbc_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      r(1:3) = r(1:3) - boxl*anint(r(1:3)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      r(:,1:3) = r(:,1:3) - boxl*anint(r(:,1:3)*boxli)
   END SUBROUTINE

   SUBROUTINE WRITE_XYZ(iu,nattached)
      USE atom_types, only: atom_name,atom
      integer,intent(in):: iu,nattached
      integer:: i
      write(iu,*) natom
      write(iu,'(a,i6)') 'Silica   nattached = ',nattached
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*Angstrom
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,it
      real(wp):: x,y,z
      character(32):: ctmp
      read(iu,*) natom
      read(iu,*) ctmp,ctmp,ctmp,it
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod

