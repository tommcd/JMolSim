MODULE pdb2xyz_conn_mod
   USE precision_mod
   USE global_vars
   USE files_mod
   USE connectivity_mod
   IMPLICIT NONE
   SAVE
   integer,parameter :: maxlnt = 132
   character(len=*),parameter :: fmtt = '(a132)'
   character(len=maxlnt) :: line
   integer :: istart,lnt
   public:: getword,getline,line,istart,lnt

CONTAINS

   SUBROUTINE pdb2xyz_conn(iu)
      integer:: iargc
      integer :: ia,ig,natom,iu,io,i,ios
      real(wp) :: rxyz(:,:)
      character(len=132) fin0,fout0,fout1
      character(len=80):: word,carg
      character(len=4):: atom_temp(:)
!
!      narg = command_argument_count()
!      write (*,*) 'number of command arguments = ', narg
!      if(narg /= 2)then
!         write(*,*)'usage :'
!         write(*,*)'  pdb2xyz_conn in.pdb out.xyz out.connect'
!      end if
!!     call get_command(ctmp, len, status)
!
!
!      call get_command_argument (1, fout0, len, status)
!      if (status /= 0) then
!         write (*,*) 'Getting command name failed with status = ', status
!         stop
!      end if
!      write (*,*) 'arg = ', ctmp(1:len)
!      read(unit=ctmp,fmt=*)ncol
!
      natom = 0
      lines: do
         call getline(iu,ios)
         if(ios /= 0) exit
         call getword(word,ios)
         select case(trim(word))
         case('HETATM')
            natom = natom + 1
            call getword(carg,ios)
            read(unit=carg,fmt=*)ia
            call getword(carg,ios)
            read(unit=carg,fmt=*)atom_temp(natom)
            call getword(carg,ios)
            read(unit=carg,fmt=*)ig
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(natom,1)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(natom,2)
            call getword(carg,ios)
            read(unit=carg,fmt=*)rxyz(natom,3)
        case('CONECT')
            call getword(carg,ios)
            read(unit=carg,fmt=*)iatom
            ii = 0
            do
               call getword(carg,ios)
               if(ios /= 0)exit
               read(unit=carg,fmt=*)jj
               if(jj == 0) cycle
               !set_proximity(iatom,0,jj)
               ii = ii + 1
               proximity(iatom,ii) = jj
            end do
        case('END')
            EXIT
        case default
            write(*,*) 'error : ',trim(word),' is not on list'
        end select
      end do lines
!
      call get_free_file_unit(io)
      open(unit=io,file= trim(fout0),status = 'unknown')
      write(io,*) natom
      write(io,*) 'Silica'
      do i = 1,natom
         write(io,'(a2,3(2x,f10.3))') trim(atom_temp(i)),x(i),y(i),z(i)
      end do
      write(io,*)
      close(io)
   END SUBROUTINE pdb2xyz_conn



   SUBROUTINE getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
   END SUBROUTINE getline

   
   SUBROUTINE getword(word,ios)
      character(len=*),intent(out) :: word
      integer,intent(out):: ios
      integer :: i,iend,i1
      ios = 0
      word(:)=''
      if(istart > lnt)then
         ios = -1
         return
      end if
      i1 = 0
      do i = istart,lnt
         if(line(i:i) /= ' ')then
            i1 = i
            exit
         end if
      end do
      if(i1 == 0) then
         ios = -2
         return
      end if
      iend = 1
      do i = i1+1,lnt
         if(line(i:i) == ' ' .or. i == lnt)then
            iend = i
            exit
         end if
      end do
      if (i1 >= lnt) iend = i1
      word(1:iend-i1+1) = line(i1:iend)
      istart = i+1
      return
   END SUBROUTINE  getword

END MODULE pdb2xyz_conn_mod

