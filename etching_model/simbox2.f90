   PROGRAM SIMBOX


      implicit none
      integer,parameter:: natom = 1536
      real*8 POSITION(natom,3),LL(3),POSITION1(natom,3)
      integer proximity(natom,27)
      real*8 CL,boxl,boxli,length
      integer x,y,z,j,i,l,m,n
      character(2) atom(natom)

      OPEN(unit=13,FILE = 'unitcell.dat',STATUS = 'OLD')
      OPEN(unit=14,FILE = 'simbox.out',STATUS = 'NEW')
      OPEN(unit=15,FILE = 'neighbor.out',STATUS = 'NEW')
!*****   READ 'input.dat'  *******************
      do  I = 1,24
         READ(13,*) POSITION(I,1),POSITION(I,2),POSITION(I,3)
      end do
      CLOSE(UNIT = 13,STATUS = 'KEEP')

!------------------------------------------------------------
      CL = 40.0d0
      boxl = 4.0d0*CL
      boxli = 1.0d0/boxl

!                     ***************************
!*********************Creating the simulation box************************
!                     ***************************
      do z = 1,4
      do y = 1,4
      do x = 1,4
         do j = 1,24
            if (j < 9) then
               write(14,*) 'Si', POSITION(j,1) + (x - 1)*CL, POSITION(j,2) + (y - 1)*CL,POSITION(j,3) + (z - 1)*CL
            else
               write(14,*) 'O', POSITION(j,1) + (x - 1)*CL, POSITION(j,2) + (y - 1)*CL, POSITION(j,3) + (z - 1)*CL
            end if
         end do
      end do
      end do
      end do

      CLOSE(UNIT = 14,STATUS = 'KEEP')
!*************************************************************************
!*************************************************************************
!*************************************************************************


!                    ****************************
!********************Making list of the neighbors*************************
!                    ****************************


      OPEN(unit=14,FILE = 'simbox.out',STATUS = 'OLD')
      do z = 1,natom
      READ(14,*) atom(z), POSITION(z,1), POSITION(z,2), POSITION(z,3)
      POSITION1(z,1) = POSITION(z,1) - 1.5d0*CL
      POSITION1(z,2) = POSITION(z,2) - 1.5d0*CL
      POSITION1(z,3) = POSITION(z,3) - 1.5d0*CL

!---------------- Periodic boundary conditions ------------------------

      POSITION(z,1) = POSITION1(z,1) - boxl*anint(POSITION1(z,1)*boxli)
      POSITION(z,2) = POSITION1(z,2) - boxl*anint(POSITION1(z,2)*boxli)
      POSITION(z,3) = POSITION1(z,3) - boxl*anint(POSITION1(z,3)*boxli)
      end do
!-----------------------------------------------------------------------
      do z = 1,natom
         i = 0
         do y = 1,natom
            do x = 1,3
               LL(x) = POSITION(y,x) - POSITION(z,x)
               LL(x) = LL(x) - boxl*anint(LL(x)*boxli)
            end do
               if (((((length(LL).LE.(CL*sqrt(2.0d0)/8.0d0)) &
               .OR.(length(LL).LE.(CL*sqrt(3.0d0/2.0d0)/4.0d0))) &
               .AND.(y.NE.z))).AND.(atom(z).NE.atom(y))) then
                i = i + 1
                proximity(z,i) = y
               end if
         end do
       end do

      do i = 1,natom
      write (15,*) atom(i),proximity(i,1),proximity(i,2),proximity(i,3),proximity(i,4)
      end do
!*************************************************************************************
!*************************************************************************************
!*************************************************************************************

stop
end

          Double precision Function length(vector1)
          DIMENSION vector1(3)
           double precision vector1
           integer y
           lenth = 0.0d0
           do  y = 1,3
             length = length + vector1(y)**2
           end do
             length = sqrt(length)
               RETURN
                 end
