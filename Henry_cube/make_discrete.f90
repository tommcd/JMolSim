PROGRAM make_discrete
      implicit none
      integer,parameter:: wp = kind(1.0d0)
      integer,parameter:: nbin = 14
      integer:: i,j,k,narg,len,status
      real(wp):: x,y,z
      real(wp):: henryk(nbin,nbin,nbin)
      integer:: ia(nbin,nbin,nbin) = 0
      character(len=132) infile,outfile,ctmp
!
      narg = command_argument_count()
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      if(narg < 2)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' infile outfile'
         stop
      end if
!
      call get_command_argument (1, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (2, outfile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', outfile(1:len)
      open(unit=15,file=trim(outfile(1:len)))

!
      do k = 1,nbin
         do j = 1,nbin
            do i = 1,nbin
               read (14,*) x,y,z,henryk(i,j,k)
               if(henryk(i,j,k) < 0.035_wp) ia(i,j,k) = 1
               write(15,'(3f9.6,i3)') x,y,z,ia(i,j,k) 
            end do
         end do
      end do
      write(15,*)
!
      do k = 1,nbin
         do j = 1,nbin               
            write(*,'(20i1)')(ia(i,j,k),i=1,nbin)
         end do
         write(*,*)
      end do
contains

   function command_argument_count()
      integer:: command_argument_count
      integer,external:: iargc
      command_argument_count = iargc()
   end function

   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
      INTEGER         , INTENT(IN)            :: NUMBER
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
      CHARACTER(LEN=1000) :: TMPVAL
      INTEGER :: IARGC
      EXTERNAL   IARGC
      IF (NUMBER < 0) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 1
          RETURN
      ELSE IF (NUMBER > IARGC()) THEN
          IF (PRESENT(VALUE )) VALUE  = ' '
          IF (PRESENT(LENGTH)) LENGTH = 0
          IF (PRESENT(STATUS)) STATUS = 2
          RETURN
      END IF
      IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
      IF (PRESENT(LENGTH)) THEN
          IF (PRESENT(VALUE)) THEN
              LENGTH = LEN_TRIM(VALUE)
          ELSE
              CALL GETARG(NUMBER,TMPVAL)
              LENGTH = LEN_TRIM(TMPVAL)
          END IF
      END IF
      IF (PRESENT(STATUS)) STATUS = 0
      RETURN
   END SUBROUTINE GET_COMMAND_ARGUMENT
end program

