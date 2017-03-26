
MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:)
   real(wp):: boxl(3),boxli(3),boxl2(3) !,boxl2i,boxl2n
   interface pbc
      module procedure  pbc_v, pbc_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
!r(1:2) = r(1:2) - boxl*anint(r(1:2)*boxli)
!r(1:2) = r(1:2) - boxl*aint(r(1:2)*boxl2i)
!if(r(1) >  boxl2) r(1) = r(1) - boxl
!if(r(1) < -boxl2) r(1) = r(1) + boxl
!if(r(2) >  boxl2) r(2) = r(2) - boxl
!if(r(2) < -boxl2) r(2) = r(2) + boxl
!if(r(1) > boxl2 ) r(1) = r(1) - boxl
!if(r(1) < boxl2n) r(1) = r(1) + boxl
!if(r(2) > boxl2 ) r(2) = r(2) - boxl
!if(r(2) < boxl2n) r(2) = r(2) + boxl
!if (rr(1) >  boxl2) then
!   rr(1) = rr(1) - boxl
!else if (rr(1) < -boxl2) then
!   rr(1) = rr(1) + boxl
!end if
!if (rr(2) >  boxl2) then
!   rr(2) = rr(2) - boxl
!else if(rr(2) < -boxl2) then
!   rr(2) = rr(2) + boxl
!end if
if(abs(r(1)) > boxl2(1)) r(1) = r(1) - sign(boxl(1),r(1))
if(abs(r(2)) > boxl2(2)) r(2) = r(2) - sign(boxl(2),r(2))
if(abs(r(3)) > boxl2(3)) r(3) = r(3) - sign(boxl(3),r(3))
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      integer:: i
!     r(:,1:2) = r(:,1:2) - boxl*anint(r(:,1:2)*boxli)
      do i=1,size(r,1)
         if(abs(r(i,1)) > boxl2(1)) r(i,1) = r(i,1) - sign(boxl(1),r(i,1))
         if(abs(r(i,2)) > boxl2(2)) r(i,2) = r(i,2) - sign(boxl(2),r(i,2))
         if(abs(r(i,3)) > boxl2(3)) r(i,3) = r(i,3) - sign(boxl(3),r(i,3))         
      end do
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
      read(iu,*) ctmp!,ctmp,ctmp,it
      it = 0
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod

