
!include 'readline_mod.f90'

module readline
    implicit none
    save
    integer,parameter :: maxlnt = 132
    character(len=*),parameter :: fmtt = '(a132)'
    character(len=maxlnt) :: line
    integer :: istart,lnt
    public:: getword,getline,line,istart,lnt
!
contains

  subroutine getline(iu,ios)
      integer,intent(in) :: iu
      integer,intent(out):: ios
      line = ''
      read(unit=iu,fmt=fmtt,iostat=ios) line
      if(ios /= 0)return
      istart = 1
      lnt = len_trim(line)
  end subroutine getline

  subroutine getword(word)
      character(len=*),intent(out) :: word
      integer :: i,iend,i1
      word(:)=' '
      if(istart > lnt)return
      do i = istart,lnt
         if(line(i:i) /= ' ')then
            i1 = i
            exit
          end if
      end do
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
  end subroutine getword
end module readline


include 'files_mod.f90'

program avg_gnu
      use readline
      use files_mod
      implicit none
      integer, parameter :: wp = selected_real_kind(15,300)
      integer:: len,status
      integer :: i,ii,j,k,narg,ncol,nfiles,iout,ia(1),ngraphs
      integer,allocatable:: iu(:),ios(:)
      character(len=132),allocatable:: files(:)
      real(wp),allocatable:: dat(:,:),sumd(:)
      character(len=132) outfile,ctmp,word
      logical:: avg_col_one = .false.
      logical:: blankline
      type graph_data
         integer:: nl,nc
         real(wp):: dat(1000,10)
      end type
      type(graph_data):: gr(10),gsum
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      if(narg < 2)then
         write(*,*)'usage :'
         write(*,*)'  avg_gnu ncols output_file file1 file2 ...'
         stop
      end if
!     call get_command(ctmp, len, status)

      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', ctmp(1:len)

      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*)ncol
      print *,'ncol = ',ncol
      call get_command_argument (2, outfile, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'arg = ',outfile(1:len)
      write (*,*),' output file = ',trim(outfile)
      call get_free_file_unit(iout)
      open(unit=iout,file=trim(outfile))

      nfiles = narg-2
      allocate(iu(nfiles),files(nfiles),ios(nfiles))
      allocate(dat(nfiles,ncol),sumd(ncol))
      i = 0
      do ii = 3, narg
         i = i + 1
         call get_command_argument (ii, files(i), len, status)
         if (status /= 0) then
            write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', ii
            stop
         end if
         write (*,*) 'arg ', ii, ' = ', files(i)(1:len)
         call get_free_file_unit(iu(i))
         open(unit=iu(i),file=trim(files(i)),status='old',iostat=ios(i))
         if (ios(i) /= 0) then
            write (*,*) 'error opening ',trim(files(i)),' ios = ',ios(i)
            stop
         end if
      end do
      write (*,*) 'command line processed'
!
      graphs: do
         ios = 0
         loop_files: do i = 1,nfiles
            gr(i)%nl = 0
            gr(i)%dat = 0.0_wp
            blankline = .false.
            lines: do
               call getline(iu(i),ios(i))  !;print *,trim(line);print *,ios(i)
               if(ios(i) /= 0) exit
               if(scan(line,'#') /= 0) cycle lines
               if(line == '')then
                  if(blankline)then
                     exit lines
                  else
                     blankline = .true.
                     cycle lines
                  end if
               else
                  blankline = .false.
                  gr(i)%nl = gr(i)%nl + 1
                  do j = 1,ncol
                     call getword(word)  !;print *,trim(word)
                     read(unit=word,fmt=*) gr(i)%dat(gr(i)%nl,j)
                  end do
               end if
            end do lines
         end do loop_files
         ngraphs = count(gr(1:nfiles)%nl > 0)
         print '(a,100i3)','gr%nl = ',gr(1:nfiles)%nl
         print *,'ngraphs = ',ngraphs
         print '(a,100i3)','ios = ',ios

         if(all(ios /= 0)) exit graphs

         gsum%nl = maxval(gr%nl)
         ia = maxloc(gr%nl)
         print '(a,i6)','maxloc(gr%nl) = ',ia
         do i = 1,nfiles
            do k = 1, gr(i)%nl
               write(*,'(10f12.6)') (gr(i)%dat(k,j),j=1,ncol)
            end do
            if(gr(i)%nl > 0)write(*,*)
         end do

         gsum%dat = 0.0
         do i = 1, nfiles
            gsum%dat = gsum%dat + gr(i)%dat
         end do
! check that the data files are consistent
         if(.not.avg_col_one)then
            do k = 1,gsum%nl
               do i = 1,nfiles
                  if(i == ia(1))cycle
                  if((gr(i)%dat(k,1) /= gr(ia(1))%dat(k,1)).and. &
                     (gr(i)%dat(k,1) /= 0.0 )) then
                     stop 'incompatible data files'
                  end if
               end do
            end do
         end if
! average the data
         gsum%dat = gsum%dat/ngraphs
         if(.not.avg_col_one)gsum%dat(:,1)=gr(ia(1))%dat(:,1)
!
         do i = 1, gsum%nl
            write(*,'(10f12.6)') (gsum%dat(i,j),j=1,ncol)
            write(iout,'(10f12.6)') (gsum%dat(i,j),j=1,ncol)
         end do
         write(*,'(/)')
         write(iout,'(/)')
       end do graphs
end program

