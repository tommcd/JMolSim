MODULE probe_mol_mod
   USE precision_mod
   implicit none
!
!  definition of probe_mol derived type
!  the allocatable component requires a standard extension to F95
!  described in the technical report ISO/IEC TR 15581
!
   TYPE probe_mol
      integer:: n
      real(wp),allocatable:: r(:,:)
      integer,allocatable:: atom(:)
      real(wp),allocatable:: rad(:),q(:)
      integer,allocatable:: nc(:),con(:,:)
   END TYPE probe_mol
!
CONTAINS

   SUBROUTINE new_probe_mol(P,n,radius,atomtype,flexible)
      type(probe_mol),intent(out):: P
      integer,intent(in):: n
      real(wp),intent(in),optional:: radius
      integer,intent(in),optional:: atomtype
      logical,intent(in),optional:: flexible
      P%n = n
      allocate(P%r(3,n),P%atom(n),P%rad(n),P%q(n))
      if (present(flexible) .and. flexible == .true.) then
         allocate(P%conect(n,4))
         P%conect = 0
      end if
      P%r = 0.0_wp
      P%atom = 0
      P%rad = 0.0_wp
      if (present(radius)) P%rad = radius
      if (present(atomtype)) P%atom = atomtype
      P%q = 0.0_wp
   END SUBROUTINE

   SUBROUTINE WRITE_PROBE_MOL(iu,P)
      USE atom_types_mod
      integer,intent(in):: iu
      type(probe_mol),intent(in):: P
      integer:: i
      write(iu,*) P%n
      write(iu,*) 'Molecule'
      do i = 1,P%n
!        write(iu,'(a6,i5,a4,2x,a3,i6,4x,3f8.3')'HETATM',i,atom_name(P%atom(i)), &
!                                               '   ',i,P%r(1:3,i)*angstrom
         write(iu,'(a2,3(1x,f14.8))') atom_name(P%atom(i)),P%r(1:3,i)*Angstrom
      end do
   END SUBROUTINE

END MODULE probe_mol_mod

