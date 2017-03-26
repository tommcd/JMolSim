
include 'precision_mod.f90'
include 'sort_mod.f90'
include 'HKNonLattice2.f90'
include 'rand_mod.f90'

PROGRAM perc2mol
      USE precision_mod
      USE HKNonLattice_mod
      USE rand_mod
      implicit none
      integer,parameter:: NX = 50
      integer,parameter:: NY = 50
      integer,parameter:: NZ = 50
      integer:: lst_top(NX*NY),ntop
      integer:: lst_span(NX*NY),nspan
      integer,parameter:: nbin = NX*NY*NZ
      integer,parameter:: full = 0
      integer,parameter:: free = 1
      real(wp):: acrit
      integer:: i,j,k,narg,length,stat,nb
      integer:: ip,im,jp,jm,kp,km,na,cc
      real(wp):: arr(NX,NY,NZ)
      integer:: abool(NX,NY,NZ)
      integer:: atnum(NX,NY,NZ)
      integer:: clnum(NX,NY,NZ)
      integer:: ncon(nbin)=0,conect(nbin,6)=0
      character(len=132):: infile,outfile,ctmp
      real:: t0,t1
!     logical:: connected
!
      call cpu_time(t0)
      call Init_HKNonLattice(nbin)
!
      narg = command_argument_count()
      call get_command_argument (0, ctmp, length, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with stat = ', stat
         stop
      end if
      if(narg < 2)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:length),' acrit  irand_seed'
         stop
      end if
!
      call get_command_argument (1, ctmp, length, stat)
      if (stat /= 0) then
         write (*,*) 'get_command_argument failed: stat = ', stat, ' arg = ', 1
         stop
      end if
      read(unit=ctmp,fmt=*) acrit
!
      call get_command_argument (2, ctmp, length, stat)
      if (stat /= 0) then
         write (*,*) 'get_command_argument failed: stat = ', stat, ' arg = ', 2
         stop
      end if
      read(unit=ctmp,fmt=*) irand_seed

!      call get_command_argument (2, ctmp, length, stat)
!      if (stat /= 0) then
!         write (*,*) 'get_command_argument failed: stat = ', stat, ' arg = ', 2
!         stop
!      end if
!      write (*,*) 'arg = ', ctmp(1:length)
!      read(unit=ctmp,fmt=*) Kend ; print *,'Kend = ',Kend
!      delK = (Kend - Kbegin)/ntry
!
!      call get_command_argument (3, infile, length, stat)
!      if (stat /= 0) then
!         write (*,*) 'get_command_argument failed: stat = ', stat, ' arg = ', 3
!         stop
!      end if
!      write (*,*) 'arg = ', infile(1:length)
!      open(unit=14,file=trim(infile(1:length)))
!
!      call get_command_argument (2, outfile, length, stat)
!      if (stat /= 0) then
!         write (*,*) 'get_command_argument failed: stat = ', stat, ' arg = ', 2
!         stop
!      end if
!      write (*,*) 'arg = ', outfile(1:length)

!
      !call random_number(arr)
      do k = 1,NZ
      do j = 1,NY
      do i = 1,NX
         arr(i,j,k) = rand()
      end do
      end do
      end do

      print *,arr(1,1,2)
!
! test code (remove later)
!      arr = 0
!      arr(2,2,:) = 1.0
!      arr(2,1,:) = 1.0
!      arr(3,2,5) = 0.5
!      arr(4,2,5) = 0.5
!      arr(4,3,5) = 0.5
!      arr(4,4,:) = 1.0
!      arr(2,10,:) = 1.0
! end test code

!
! Where array arr > acrit the abool array is labeled 'free'
! and an 'atom' na is placed in the atnum array.
!
      na = 0
      do k = 1,NZ
      do j = 1,NY
      do i = 1,NX
         if (arr(i,j,k) > acrit) then
            na = na + 1
            abool(i,j,k) = free
            atnum(i,j,k) = na
         else
            abool(i,j,k) = full
            atnum(i,j,k) = 0
         end if
      end do
      end do
      end do
!
      print *,'acrit = ',acrit
      write(*,*) 'bins with a > acrit ',count(abool == free)
!
! Loop over the free elements and connect the 'atoms'
! the array is periodic in XY directions but NOT Z
!
      conect = 0
      na = 0
      do k = 1,NZ
         kp = k + 1
         km = k - 1
         !if (kp > NZ) kp = 1
         !if (km < 1) km = NZ
         do j = 1,NY
            jp = j + 1
            if (jp > NY) jp = 1
            jm = j - 1
            if (jm < 1) jm = NY
            do i = 1,NX
               if (abool(i,j,k) == full) CYCLE
               na = na + 1
               ip = i + 1
               if (ip > NX) ip = 1
               im = i - 1
               if (im < 1) im = NX
               ! connectivity
               nb = 0
               if (atnum(ip,j,k) /= 0) then
                  nb = nb + 1
                  conect(na,nb) = atnum(ip,j,k)
               end if
               if (atnum(im,j,k) /= 0) then
                  nb = nb + 1
                  conect(na,nb) = atnum(im,j,k)
               end if
               if (atnum(i,jp,k) /= 0) then
                  nb = nb + 1
                  conect(na,nb) = atnum(i,jp,k)
               end if
               if (atnum(i,jm,k) /= 0) then
                  nb = nb + 1
                  conect(na,nb) = atnum(i,jm,k)
               end if
               if (kp <= NZ) then
                  if (atnum(i,j,kp) /= 0) then
                     nb = nb + 1
                     conect(na,nb) = atnum(i,j,kp)
                  end if
               end if
               if (km >= 1) then
                  if (atnum(i,j,km) /= 0) then
                     nb = nb + 1
                     conect(na,nb) = atnum(i,j,km)
                  end if
               end if
               ncon(na) = nb
            end do         
         end do
      end do
!
! Now use the Non-Lattice version of the Hochen-Kopelman cluster
! algorithm to label the clusters of connectd 'atoms' i.e. 
! connected 'free' elements of the array abool. atomL(i)
! stores the label (cluster number) of atom i
!
      if (na /= 0) then
         call HKNonLattice(na,conect,n_cluster,atomL)
!        call analyse_cluster(na)
!        call print_cluster(na,6)
      else
         stop 'no nodes'
      end if

!
! If abool(i,j,k) is 'free' (& so contains an 'atom') then
! clnum(i,j,k) is assigned it's cluster number. Otherwise
! it is set to zero.
!
      ! clnum = 0
      ! WHERE(atnum /= 0) clnum = atomL(atnum)
      do k = 1,NZ
      do j = 1,NY
      do i = 1,NX
         if (atnum(i,j,k) == 0)then
            clnum(i,j,k) = 0
            CYCLE
         else
            clnum(i,j,k) = atomL(atnum(i,j,k))
         end if
      end do         
      end do
      end do

!
! Check each element of the top z-slice of the array i.e. clnum(:,:,NZ)
! and store the ntop unique clusters in a list (lst_top).
!
      ntop = 0
      lst_top = 0
      do j = 1,NY
      do i = 1,NX
         cc = clnum(i,j,NZ)
         if (cc == 0) cycle
         if (.not.in_list(lst_top,ntop,cc)) then
            ntop = ntop + 1
            lst_top(ntop) = cc
         end if
      end do         
      end do

!     print *,'ntop = ',ntop
!     print *,'lst_top: ',lst_top(1:ntop)
!
! Check each element of the bottom z-slice of the array i.e. clnum(:,:,1)
! and store the nspan unique clusters which are also in the top layer
! (i.e. the z-spanning clusters) in a list (lst_span).
!
      nspan = 0
      lst_span = 0
      do j = 1,NY
      do i = 1,NX
         cc = clnum(i,j,1)
         if (cc == 0) cycle
         if (in_list(lst_top,ntop,cc)) then
            if (.not.in_list(lst_span,nspan,cc)) then
               nspan = nspan + 1
               lst_span(nspan) = cc
            end if
         end if
      end do         
      end do

      print *,' there are ',nspan,' z-spanning clusters:'
      do i = 1,nspan
         print *,lst_span(i)
      end do
      call cpu_time(t1)
      print *,'time taken ',t1-t0, ' s'
      
contains

   pure function in_list(Ls,n,cc)
      logical:: in_list
      integer,intent(in):: Ls(:),n,cc
      integer:: k
      in_list = .false.
      do k = 1,n
         if (cc == Ls(k)) then
            in_list = .true.
            return
         end if
      end do
   end function

!   function command_argument_count()
!      integer:: command_argument_count
!      integer,external:: iargc
!      command_argument_count = iargc()
!   end function
!
!   SUBROUTINE GET_COMMAND_ARGUMENT(NUMBER,VALUE,LENGTH,STATUS)
!      INTEGER         , INTENT(IN)            :: NUMBER
!      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VALUE
!      INTEGER         , INTENT(OUT), OPTIONAL :: LENGTH
!      INTEGER         , INTENT(OUT), OPTIONAL :: STATUS
!      CHARACTER(LEN=1000) :: TMPVAL
!      INTEGER :: IARGC
!      EXTERNAL   IARGC
!      IF (NUMBER < 0) THEN
!          IF (PRESENT(VALUE )) VALUE  = ' '
!          IF (PRESENT(LENGTH)) LENGTH = 0
!          IF (PRESENT(STATUS)) STATUS = 1
!          RETURN
!      ELSE IF (NUMBER > IARGC()) THEN
!          IF (PRESENT(VALUE )) VALUE  = ' '
!          IF (PRESENT(LENGTH)) LENGTH = 0
!          IF (PRESENT(STATUS)) STATUS = 2
!          RETURN
!      END IF
!      IF (PRESENT(VALUE)) CALL GETARG(NUMBER,VALUE)
!      IF (PRESENT(LENGTH)) THEN
!          IF (PRESENT(VALUE)) THEN
!              LENGTH = LEN_TRIM(VALUE)
!          ELSE
!              CALL GETARG(NUMBER,TMPVAL)
!              LENGTH = LEN_TRIM(TMPVAL)
!          END IF
!      END IF
!      IF (PRESENT(STATUS)) STATUS = 0
!      RETURN
!   END SUBROUTINE GET_COMMAND_ARGUMENT
end program

