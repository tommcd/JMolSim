SUBROUTINE Energy_Minimization_DPD
      USE CommonBlocks
      USE CommonDPD
      USE Eminparam
      implicit none
!
! minimize the energy
! Steepest Decent Algorithm
!
      real(8), parameter :: dRmaxDP = 0.1
      real(8), dimension(3,N) :: Rtmp
      real(8) :: sumF2, RmsF, RmsF_o
      real(8) :: deltaE, Energy, Energy_Pre
      real(8) :: deltaR
      integer :: Itry, i

      print *, 'Enegy minimization has just started'

      deltaR = dRmaxDP

      call ForceDPD_Minim_Original

      Energy_Pre = PotDP
      write( 6,*) 'Initial Energy             : ',Energy_Pre
      write(11,*) 'Initial Energy             : ',Energy_Pre

      sumF2 = 0.d0
      do i = 1 , N
         sumF2 = sumF2 + dot_product( FrcDP(:,i), FrcDP(:,i) )
      end do
      RmsF = sqrt( sumF2 / Nf )

!
! Calculate Maximum Downhill Gradient
!
      do Itry = 1 , MinTry

         RmsF_o = RmsF
         Rtmp   = R
         FrcDPt = FrcDP

         R = R + deltaR * FrcDP / RmsF

         do i = 1 , N
            R(:,i) = R(:,i) - nint( R(:,i) * InvCL ) * CellL
         end do

         call ForceDPD_Minim_Original

         Energy = PotDP

         deltaE = Energy - Energy_Pre

         sumF2 = 0.d0
         do i = 1 , N
            sumF2 = sumF2 + dot_product( FrcDP(:,i), FrcDP(:,i) )
         end do
         RmsF = sqrt( sumF2 / Nf )

         write(*,'(a,i3,a,e12.4)') 'Step = ',Itry,' :  Energy = ', Energy

         if ( ( abs(deltaE/Energy) < dev_relative ).or.(Energy == 0.) ) exit

         if ( Energy < Energy_Pre ) then
            Energy_Pre = Energy
            deltaR     = deltaR*1.2d0
            deltaR     = min(deltaR , dRmaxDP)
         else
            FrcDP  = FrcDPt
            R      = Rtmp
            RmsF   = RmsF_o
            deltaR = deltaR*0.5d0
         end if

      end do

      write(11,*) 'Final Energy               : ',Energy

      print *, 'Enegy minimization has just finished'
      print *, 'Final Energy               : ',Energy

END SUBROUTINE Energy_Minimization_DPD

