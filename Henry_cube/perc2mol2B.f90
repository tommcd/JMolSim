
include 'sort_mod.f90'
include 'HKNonLattice2.f90'

PROGRAM perc2mol
      USE HKNonLattice_mod
      implicit none
      integer,parameter:: wp = kind(1.0d0)
      integer,parameter:: N = 14
      integer:: lst_top(N*N),ntop
      integer:: lst_span(N*N),nspan
      integer,parameter:: nbin = N**3
      real(wp),parameter:: angstrom = 10.0_wp
      integer,parameter:: full = 0,free = 1
      real(wp):: Kcrit = 0.035_wp
      integer:: i,j,k,narg,len,status,nb,iu
      integer:: ip,im,jp,jm,kp,km,na,cc
      real(wp):: x,y,z
      real(wp):: henryk(N,N,N)
      real(wp):: rxyz(nbin,3),hk(nbin)
      integer:: a(N,N,N) = full
      integer:: b(N,N,N) = 0,c(N,N,N) = 0
      integer:: nc(nbin)=0,proximity(nbin,6)=0
      character(len=132):: infile,outfile,ctmp
      real(wp):: boxl,boxl2 !,boxli
!     logical:: connected
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
!     boxli = 1.0_wp/boxl
      call Init_HKNonLattice(nbin)
!
      narg = command_argument_count()
      call get_command_argument (0, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      if(narg < 3)then
         write(*,*)'usage :'
         write(*,*) ctmp(1:len),' Kcrit infile outfile'
         stop
      end if
!
      call get_command_argument (1, ctmp, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', ctmp(1:len)
      read(unit=ctmp,fmt=*) Kcrit ; print *,'Kcrit = ',Kcrit
!
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, outfile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', outfile(1:len)
      iu = 15
      open(unit=iu,file=trim(outfile(1:len)))

!
      na = 0
      do k = 1,N
         do j = 1,N
            do i = 1,N
               read (14,*) x,y,z,henryk(i,j,k)
!print *,x,y,z,henryk(i,j,k)
               if (henryk(i,j,k) > Kcrit) then
                  na = na + 1
                  a(i,j,k) = free
                  rxyz(na,:) = (/ x,y,z /) - (/ boxl2,boxl2,boxl2 /)
                  hk(na) = henryk(i,j,k)/Kcrit
               end if
            end do
         end do
      end do
!
      do k = 1,N
         do j = 1,N               
            write(*,'(20i1)')(a(i,j,k),i=1,N)
         end do
         write(*,*)
      end do
      write(*,*) 'bins with K > Kcrit ',count(a == free)

      b = 0
      na = 0
      do k = 1,N
      do j = 1,N
      do i = 1,N
         ! ib = i + (j-1)*N + (k-1)*N*N
         if (a(i,j,k) == full) CYCLE
         na = na + 1
         b(i,j,k) = na
      end do         
      end do
      end do
!
      proximity = 0
      na = 0
      do k = 1,N
         kp = k + 1
         if (kp > N) kp = 1
         km = k - 1
         if (km < 1) km = N
         do j = 1,N
            jp = j + 1
            if (jp > N) jp = 1
            jm = j - 1
            if (jm < 1) jm = N
            do i = 1,N
               if (a(i,j,k) == full) CYCLE
               na = na + 1
               ip = i + 1
               if (ip > N) ip = 1
               im = i - 1
               if (im < 1) im = N
               ! connectivity
               nb = 0
               if (b(ip,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(ip,j,k)
               end if
               if (b(im,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(im,j,k)
               end if
               if (b(i,jp,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jp,k)
               end if
               if (b(i,jm,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jm,k)
               end if
               if (b(i,j,kp) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,j,kp)
               end if
               if (b(i,j,km) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,j,km)
               end if
               nc(na) = nb
            end do         
         end do
      end do
!
      if (na /= 0) then
         call HKNonLattice(na,proximity,n_cluster,atomL)
         call analyse_cluster(na)
         call print_cluster(na,6)
      end if
!
! nc(i) is the number of elements in cluster i
! the array cluster(i,j) contains the j'th element of cluster i

      do ic = 1, n_cluster
         ! find the element which is closest to the centre of cluster i
         rtmp = 0.0_wp
         rr = 0.0_wp
         do ia = 1,nc(ic)
            j = cluster(ic,ia)
            rtmp = rtmp + rxyz(j,:)
            rr(ia,:) = rxyz(j,:)
         end do
         rmean = rtmp/nc(ic)
         rtmp = 0.0_wp
         do ia = 1,nc(ic)
            j = cluster(ic,ia)
            rtmp = rtmp + rxyz(j,:)
         end do
      end do

! remove z periodicity
      proximity = 0
      na = 0
      do k = 1,N
         kp = k + 1
         km = k - 1
         !if (kp > N) kp = 1
         !if (km < 1) km = N
         do j = 1,N
            jp = j + 1
            if (jp > N) jp = 1
            jm = j - 1
            if (jm < 1) jm = N
            do i = 1,N
               if (a(i,j,k) == full) CYCLE
               na = na + 1
               ip = i + 1
               if (ip > N) ip = 1
               im = i - 1
               if (im < 1) im = N
               ! connectivity
               nb = 0
               if (b(ip,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(ip,j,k)
               end if
               if (b(im,j,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(im,j,k)
               end if
               if (b(i,jp,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jp,k)
               end if
               if (b(i,jm,k) /= 0) then
                  nb = nb + 1
                  proximity(na,nb) = b(i,jm,k)
               end if
               if (kp <= N) then
                  if (b(i,j,kp) /= 0) then
                     nb = nb + 1
                     proximity(na,nb) = b(i,j,kp)
                  end if
               end if
               if (km >= 1) then
                  if (b(i,j,km) /= 0) then
                     nb = nb + 1
                     proximity(na,nb) = b(i,j,km)
                  end if
               end if
               nc(na) = nb
            end do         
         end do
      end do
!
      print *,'removed z perodicity'
      if (na /= 0) then
         call HKNonLattice(na,proximity,n_cluster,atomL)
         call analyse_cluster(na)
         call print_cluster(na,6)
      else
         stop 'no nodes'
      end if

      c = 0
      do k = 1,N
      do j = 1,N
      do i = 1,N
         if (b(i,j,k) == 0) CYCLE
         c(i,j,k) = atomL(b(i,j,k))
      end do         
      end do
      end do

      ntop = 0
      lst_top = 0
      do j = 1,N
      do i = 1,N
         cc = c(i,j,1)
         if (cc == 0) cycle
         if (.not.in_list(lst_top,ntop,cc)) then
            ntop = ntop + 1
            lst_top(ntop) = cc
         end if
      end do         
      end do

      print *,'ntop = ',ntop
      print *,'lst_top: ',lst_top(1:ntop)

      nspan = 0
      lst_span = 0
      do j = 1,N
      do i = 1,N
         cc = c(i,j,N)
         if (cc == 0) cycle
         if (in_list(lst_top,ntop,cc)) then
            if (.not.in_list(lst_span,nspan,cc)) then
               nspan = nspan + 1
               lst_span(nspan) = cc
            end if
         end if
      end do         
      end do

      print *,'there are ',nspan,' spanning clusters in the z-direction:'
      do i = 1,nspan
         print *,lst_span(i)
      end do

!      outer_loop: do k = 2,N
!         connected = .false.
!         do j = 1,N
!         do i = 1,N
!            cc = c(i,j,k)
!            if (cc == 0) cycle
!            if (in_list(lst_top,ntop,cc)) then
!               connected = .true.
!               cycle outer_loop
!            end if
!         end do         
!         end do
!         if (.not.connected) exit outer_loop
!      end do outer_loop


!
      do i = 1,na
         write(iu,110)'HETATM',i,' C  ','   ',i,rxyz(i,1:3)*angstrom,hk(i),hk(i)
      end do
      do i = 1,na
         if ( all(proximity(i,:) == 0) ) cycle
         write(iu,"(a6,7i5)")'CONECT',i,(proximity(i,j),j = 1,nc(i))
      end do
      write(iu,110)'HETATM',99991,'AR','   ',99991,-boxl2*angstrom,-boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99992,'AR','   ',99992,-boxl2*angstrom, boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2*angstrom,-boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2*angstrom, boxl2*angstrom,-boxl2*angstrom
      write(iu,110)'HETATM',99995,'AR','   ',99995,-boxl2*angstrom,-boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99996,'AR','   ',99996,-boxl2*angstrom, boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2*angstrom,-boxl2*angstrom, boxl2*angstrom
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2*angstrom, boxl2*angstrom, boxl2*angstrom
      write(iu,"(a6,5i5)")'CONECT',99991,99992
      write(iu,"(a6,5i5)")'CONECT',99991,99993
      write(iu,"(a6,5i5)")'CONECT',99991,99995
      write(iu,"(a6,5i5)")'CONECT',99992,99994
      write(iu,"(a6,5i5)")'CONECT',99992,99996
      write(iu,"(a6,5i5)")'CONECT',99993,99994
      write(iu,"(a6,5i5)")'CONECT',99993,99997
      write(iu,"(a6,5i5)")'CONECT',99994,99998
      write(iu,"(a6,5i5)")'CONECT',99995,99996
      write(iu,"(a6,5i5)")'CONECT',99995,99997
      write(iu,"(a6,5i5)")'CONECT',99996,99998
      write(iu,"(a6,5i5)")'CONECT',99997,99998
      write(iu,'(a/)')'END'
      close(iu,status='keep')
110   format(a6,i5,a4,2x,a3,i6,4x,3f8.3,2f6.2)
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

