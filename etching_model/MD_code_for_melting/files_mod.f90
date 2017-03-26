
MODULE files_mod
   public get_free_file_unit, myflush
   contains

   SUBROUTINE get_free_file_unit(iu)
      implicit none
      character(len=*),parameter :: sub_name='get_free_file_unit'
      integer,intent(out):: iu
      integer,parameter:: istart=11
      integer,parameter:: unit_max=100000
      integer:: ios
      logical:: lopen,lexists
      do iu = istart,unit_max
         inquire(unit=iu, iostat=ios, exist=lexists, opened=lopen)
         if (ios == 0) then
            if (lexists .and. .not.lopen) return
         else
            write(*,*) sub_name,': error, iostat = ',ios
            stop
         end if
      end do
      write(*,*) sub_name,': error, no free units'
      stop
   END SUBROUTINE get_free_file_unit

   SUBROUTINE myflush(iu)
      integer,intent(in):: iu
      character(132):: fn
      inquire(unit=iu,name=fn)
      close(iu)
      open(unit=iu,file=trim(fn),position='append')
   END SUBROUTINE myflush

END MODULE files_mod

