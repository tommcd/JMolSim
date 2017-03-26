   SUBROUTINE Steepest_Descent(nstep,dRmax,nat,Rxyz,Fxyz)
      integer,intent(in):: nstep
      real(wp),intent(in):: dRmax
      integer,intent(in):: nat
      real(wp),intent(inout):: Rxyz(:,:),Fxyz(:,:)
      ! real(wp), parameter :: dRmax = 0.1_wp
!
! Minimize the energy
! Steepest Descent Algorithm
!
      real(wp):: Rxyzo(nat,3),Fxyzo(nat,3)
!     real(wp),dimension(size(rxyz,1),size(rxyz,2)):: Rxyzo,Fxyzo
      real(wp):: Unew, Uold
      real(wp):: deltaR, delRF
      integer :: itry

      deltaR = dRmax
      call Force(nat,Rxyz,Fxyz,Uold)
      write (*,*) 'Initial Energy        : ',Uold

      do itry = 1, nstep

         Rxyzo = Rxyz
         Fxyzo = Fxyz

         delRF = deltaR/maxval(abs(Fxyz))

         Rxyz = Rxyz + delRF*Fxyz
         call pbc(Rxyz)

         call Force(nat,Rxyz,Fxyz,Unew)

         write(*,'(a,i3,a,e12.4)') 'Step = ',Itry,' :  Energy = ', Unew

         if (Unew < Uold) then
            Uold = Unew
            deltaR = deltaR*1.2_wp
            deltaR = min(deltaR , dRmax)
         else
            Fxyz = Fxyzo
            Rxyz = Rxyzo
            deltaR = deltaR*0.5_wp
         end if

      end do

      write (*,*) 'Final Energy          : ',Uold
   END SUBROUTINE Steepest_Descent

