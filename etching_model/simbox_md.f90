         PROGRAM SILICA
C
C
C       This program simulates the silica substance
C      using (BMH) pair potential with adding three - body
C      interections by the Stillinger and Weber potential
C
C

C
C*****   READ 'input.dat'  *******************

      implicit none
      integer,parameter:: natomt = 1536
      real(wp) POSITION(1536,3),LL(3)
      integer proximity(1536,27)
      real(wp) CL,POSITION,boxl,boxli,LL ,length
      integer x,y,z,j,i,proximity
      character(2) atom(1536)

      OPEN(unit=13,FILE = 'unitcell.dat',STATUS = 'OLD')
      OPEN(unit=14,FILE = 'simbox.out',STATUS = 'NEW')
      OPEN(unit=15,FILE = 'neighbor.out',STATUS = 'NEW')

          do 12 I = 1,24
      READ(13,*) POSITION(I,1),POSITION(I,2),POSITION(I,3)
12         continue
      CLOSE(UNIT = 13,STATUS = 'KEEP')

C - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - --
      CL = 40.0d0
      boxl = 4.0d0*CL
      boxli = 1.0d0/boxl

C                     ***************************
C*********************Creating the simulation box************************
C                     ***************************
       do z = 1,4
          do y = 1,4
             do x = 1,4
                   do j = 1,24
          if (j.LT.9) then
        write(14,*) 'Si', POSITION(j,1) + (x - 1)*CL, POSITION(j,2) + (y - 1)*CL,
     &        POSITION(j,3) + (z - 1)*CL
          else
        write(14,*) 'O', POSITION(j,1) + (x - 1)*CL, POSITION(j,2) + (y - 1)*CL,
     &        POSITION(j,3) + (z - 1)*CL
          end if

                   end do
             end do
          end do
       end do

         CLOSE(UNIT = 14,STATUS = 'KEEP')
C*************************************************************************
C*************************************************************************
C*************************************************************************


C                    ****************************
C********************Making list of the neighbors*************************
C                    ****************************


      OPEN(unit=14,FILE = 'simbox.out',STATUS = 'OLD')

              do z = 1,1536

        READ(14,*) atom(z), POSITION(z,1), POSITION(z,2), POSITION(z,3)

        POSITION(z,1) = POSITION(z,1) - 1.5d0*CL
        POSITION(z,2) = POSITION(z,2) - 1.5d0*CL
        POSITION(z,3) = POSITION(z,3) - 1.5d0*CL

C - -- - -- - -- - -- - --- Periodic boundary conditions - - -- - -- - -- - -- - -- - -- - -- - -

       POSITION(z,1) = POSITION(z,1) - boxl*anint(POSITION(z,1)*boxli)
       POSITION(z,2) = POSITION(z,2) - boxl*anint(POSITION(z,2)*boxli)
       POSITION(z,3) = POSITION(z,3) - boxl*anint(POSITION(z,3)*boxli)

C - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -
              end do

         do z = 1,1536
               i = 0
             do y = 1,1536
               do x = 1,3
                 LL(x) = POSITION(y,x) - POSITION(z,x)
                 LL(x) = LL(x) - boxl*anint(LL(x)*boxli)
               end do
          if (((((length(LL).LE.(CL*sqrt(2.0d0)/8.0d0))
     & .OR.(length(LL).LE.(CL*sqrt(3.0d0/2.0d0)/4.0d0)))
     & .AND.(y.NE.z))).AND.(atom(z).NE.atom(y))) then
                i = i + 1
                proximity(z,i) = y
            end if
               end do
         end do

          do i = 1,1536
            write (15,*) atom(i),proximity(i,1),proximity(i,2),
     &                proximity(i,3),proximity(i,4)
          end do

C*************************************************************************
C*************************************************************************
C*************************************************************************


C                    ****************************
C********************   Performance of the MD    *************************
C                    ****************************

C                              TIME OF CALCULATION
 50              if ((time).LE.(600.0d0*(1.0E-12))) then
        write(*,*) time
      number5 = number5 + 1
C                     USING THE VERLET - POSITION ALGORITHM FOR OXIGEN ATOM

                        do 20 I = 1,1024
                             do 21 J = 1,3
      POSITIONNEW(I,J) = POSITION(I,J) + DT*vel(i,j) + (ss)*DT**2*FR(I,J)/mOx

      POSITIONNEW(i,j) = POSITIONNEW(i,j)-
     & boxl*DINT(POSITIONNEW(i,j)*boxli)

        POSITION(I,J) = POSITIONNEW(I,J)
  21                         continue
  20                    continue

C                     USING THE VERLET - POSITION ALGORITHM FOR SILICON ATOM

                        do 22 I = 1025,1536
                             do 23 J = 1,3
      POSITIONNEW(I,J) = POSITION(I,J) + DT*vel(i,j) + (ss)*DT**2*FR(I,J)/mSi

      POSITIONNEW(i,j) = POSITIONNEW(i,j)-
     & boxl*DINT(POSITIONNEW(i,j)*boxli)

        POSITION(I,J) = POSITIONNEW(I,J)

  23                         continue
  22                    continue

        call FORCE
C                     USING THE VERLET - VELOCITY ALGORITHM FOR OXIGEN ATOM

                        do 40 I = 1,1024
                             do 41 J = 1,3
      VELNEW(I,J) = VEL(I,J) + (0.5d0)*DT*(FROLD(i,j)/mOx + FR(i,j)/mOx)
        VEL(I,J) = VELNEW(I,J)
  41                            continue
  40                       continue

C                     USING THE VERLET - VELOCITY ALGORITHM FOR SILICON ATOM

                        do 46 I = 1025,1536
                             do 47 J = 1,3
      VELNEW(I,J) = VEL(I,J) + (0.5d0)*DT*(FROLD(i,j)/mSi + FR(i,j)/mSi)
        VEL(I,J) = VELNEW(I,J)
  47                            continue
  46                       continue

      do 42 I = 1,1536
      do 43 J = 1,3
      FROLD(i,j) = FR(i,j)
 43    continue
 42    continue


C*************************************************************************
C*************************************************************************
C*************************************************************************


            CLOSE(UNIT = 15,STATUS = 'KEEP')
              CLOSE(UNIT = 14,STATUS = 'KEEP')

       STOP
       END

C****************** ! SUBPROGRAMS ! ****************************************
C - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - --
C - -- - -- - -- - -- Function subprogram for the absolut value of the vector
C - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - -- - --
          Double precision Function length(vector1)
      implicit none
          real(8) vector1(3)
           integer y
           length = 0.0d0
           do 10 y = 1,3
             length = length + vector1(y)**2
 10           continue
             length = sqrt(length)
               RETURN
                 end
