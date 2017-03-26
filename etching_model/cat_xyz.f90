
include 'precision_mod.f90'
include 'command_line_mod.f90'
include 'command_arg_mod.f90'
include 'files_mod.f90'

program cat_xyz
      USE precision_mod
      USE command_arg_mod
      USE files_mod
      implicit none
      integer:: iu,io,narg,iframe,nframe,natom
      real(wp),allocatable:: rxyz(:,:)
      character(len=4),allocatable:: atom_name(:)
      character(len=132):: base_name,outfile,ctmp
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call get_com_arg(0, ctmp, "command")
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) trim(ctmp),' natom nframe base_name outfile'
         write(*,*) 'natom = no. of atoms'
         write(*,*) 'nframe = no. of frames'
         write(*,*) 'base_name, frame files are named basename_ddddd.xyz'
         write(*,*) 'outfile = outputfile'
         stop 'error: incorrect no. of command line arguments'
      end if
!
      call get_com_arg(1, natom, "natom")
      call get_com_arg(2, nframe, "nframe")
      call get_com_arg(3, base_name, "base_name")
      call get_com_arg(4, outfile, "outfile")
!
      call get_free_file_unit(io)
      open(unit=io,file=trim(outfile))
      allocate(rxyz(natom,3),atom_name(natom))
! write frame zero
      call open_frame_file(iu,trim(base_name),0)
      call read_xyz(iu,natom,rxyz,atom_name)
      write(io,*) natom
      write(io,*)'xyz movie'
      call write_xyz_frame(io,natom,rxyz,atom_name)
!
      main_loop:do iframe = 1,nframe
         call open_frame_file(iu,trim(base_name),iframe)
         call read_xyz(iu,natom,rxyz,atom_name)
         call write_xyz_frame(io,natom,rxyz,atom_name)
         close(iu)
      end do main_loop
      close(io)
CONTAINS

   SUBROUTINE read_xyz(iu,natom,rxyz,atom_name)
      integer,intent(in):: iu,natom
      real(wp),intent(out):: rxyz(:,:)
      character(len=*),intent(out):: atom_name(:)
      real(wp):: rtmp,x,y,z
      integer:: i,itmp
      character(32):: ctmp
      read(iu,*) itmp
      read(iu,*) ctmp,rtmp
      do i = 1,natom
ctmp(1:2) = "  "
100      read(iu,*) ctmp,x,y,z
if (ctmp(1:2) == 'Si'.or.ctmp(1:2) == 'O '.or.ctmp(1:2) == 'OH') GOTO 100
         rxyz(i,:) = (/ x,y,z /)
         atom_name(i) = trim(ctmp)
      end do
   END SUBROUTINE read_xyz

   SUBROUTINE write_xyz_frame(iu,natom,rxyz,atom_name)
      integer,intent(in):: iu,natom
      real(wp),intent(in):: rxyz(:,:)
      character(len=*),intent(in):: atom_name(:)
      integer:: i,j
      do i = 1,natom
         write(iu,'(a,3(1x,f12.6))') trim(atom_name(i)),(rxyz(i,j),j = 1,3)
      end do
      write(iu,'(/)')
   END SUBROUTINE write_xyz_frame

   SUBROUTINE open_frame_file(iu,base_name,iframe)
      USE files_mod
      integer,intent(in):: iframe
      integer,intent(out):: iu
      character(len=*),intent(in):: base_name
      character(len=32):: ctmp
      character(len=132):: frame_file
      integer:: ios
      logical:: lexists,lopen
      call get_free_file_unit(iu)
      write(unit=ctmp,fmt='(i5.5)') iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.xyz'
      inquire(file=trim(frame_file),iostat=ios,exist = lexists,opened = lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit=iu,file=trim(adjustl(frame_file)),form = 'formatted')
   END SUBROUTINE open_frame_file

end program cat_xyz

