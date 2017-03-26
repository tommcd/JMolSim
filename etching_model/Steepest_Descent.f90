
   SUBROUTINE Steepest_Descent1(nstep,dRmax,nat,Rxyz,Fxyz)
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
   END SUBROUTINE Steepest_Descent1


   SUBROUTINE Steepest_Descent2(nstep,dRmax,nat,Rxyz,Fxyz)
      integer,intent(in):: nstep
      real(wp),intent(in):: dRmax
      integer,intent(in):: nat
      real(wp),intent(inout):: Rxyz(:,:),Fxyz(:,:)
      real(wp), parameter :: dev_rel = 1.0e-4_wp
      ! real(wp), parameter :: dRmax = 0.1_wp
!
! Minimize the energy
! Steepest Descent Algorithm
!
      real(wp):: Rxyzo(nat,3),Fxyzo(nat,3)
!     real(wp),dimension(size(rxyz,1),size(rxyz,2)):: Rxyzo,Fxyzo
      real(wp):: Unew, Uold
      real(wp):: deltaR, RmsF, RmsF_o
      integer :: itry

      deltaR = dRmax
      call Force(nat,Rxyz,Fxyz,Uold)
      RmsF = rms(Fxyz,N)
      write (*,*) 'Initial Energy        : ',Uold

      do itry = 1, nstep

         RmsF_o = RmsF
         Rxyzo = Rxyz
         Fxyzo = Fxyz

         delRF = deltaR/RmsF

         Rxyz = Rxyz + delRF*Fxyz
         call pbc(Rxyz)

         call Force(nat,Rxyz,Fxyz,Unew)

         RmsF = rms(Fxyz,N)

         write(*,'(a,i3,a,e12.4)') 'Step = ',Itry,' :  Energy = ', Unew

         if ( (abs((Unew - Uold)/Unew) < dev_rel) ) exit

         if (Unew < Uold) then
            Uold = Unew
            deltaR = deltaR*1.2_wp
            deltaR = min(deltaR , dRmax)
         else
            Fxyz = Fxyzo
            Rxyz = Rxyzo
            RmsF = RmsF_o
            deltaR = deltaR*0.5_wp
         end if

      end do

      write (*,*) 'Final Energy          : ',Unew
   CONTAINS

   pure function rms(f,n)
      real(wp):: rms
      integer,intent(in):: n
      real(wp),intent(in):: f(:,:)
      real(wp):: sum2
      integer:: i
 !    rms = sqrt(sum( Fxyz*Fxyz ) )
 !    rms = sqrt(sum( Fxyz(1:n,1:3)**2 ) )
      sum2 = 0.0_wp
      do i = 1 , n
         sum2 = sum2 + dot_product( Fxyz(i,:), Fxyz(i,:) )
      end do
      rms = sqrt(sum2/n)
   end function rms

   END SUBROUTINE Steepest_Descent2


