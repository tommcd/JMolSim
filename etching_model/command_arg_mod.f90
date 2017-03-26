
MODULE command_arg_mod
   USE precision_mod
   USE command_line_mod
   implicit none
   !
   ! Gets a command argument and assigns it to a variable.
   ! It checks the status of the argument and writes a record to stdout.
   !
   interface get_com_arg
      module procedure get_com_arg_int, get_com_arg_real, &
                       get_com_arg_logical,get_com_arg_string
   end interface
CONTAINS

   SUBROUTINE get_com_arg_int(ia,var,var_name)
      integer,intent(in):: ia
      integer,intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_real(ia,var,var_name)
      integer,intent(in):: ia
      real(wp),intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_logical(ia,var,var_name)
      integer,intent(in):: ia
      logical,intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE

   SUBROUTINE get_com_arg_string(ia,var,var_name)
      integer,intent(in):: ia
      character(*),intent(out):: var
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc,stat
      call get_command_argument(ia, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', ia
         stop
      end if
!      write (*,*) 'arg = ', carg(1:nc)
      var = ''
      var(1:nc) = carg(1:nc)
      write (*,*) var_name,' = ',trim(var)
   END SUBROUTINE

END MODULE command_arg_mod

