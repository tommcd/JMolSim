
MODULE RDF_MOD
   USE precision_mod, only: wp
   USE global_vars, only: natom,nseed
   USE coordinates_mod
   implicit none
   integer,parameter,private:: nbin=200
   real(wp),private:: RDF(nbin)
   real(wp),private:: DL1
   PUBLIC:: RDF_CALC, RDF_PRINT
CONTAINS

   SUBROUTINE RDF_CALC
      real(wp):: rv(3),r,r2,rmax,rmax2
      integer:: i,j,k
      rmax = boxl*0.5_wp
      rmax2 = rmax*rmax
      DL1 = rmax/real(nbin,wp)
      RDF = 0.0_wp
      do i = nseed+1,natom-1
         do j = i + 1,natom
            rv(1:3) = rxyz(j,1:3) - rxyz(i,1:3)
            call pbc(rv)
            r2 = dot_product(rv,rv)
            if (r2 >= rmax2) cycle
            r = sqrt(r2)
            if (r < 0.01_wp) then
               write( *,'(a,f0.6,2i6)') '# RDF_CALC: ',r,i,j
            end if
            k = int(r/DL1)+1
            RDF(k) = RDF(k) + 1.0_wp
         end do
      end do
   END SUBROUTINE RDF_CALC

   SUBROUTINE RDF_PRINT(io,nattached)
      integer,intent(in):: io,nattached
      integer:: i
      write(io,*)'# nattached = ',nattached
      write(io,*)'# (natom-nseed) = ',natom-nseed
      do i = 1,nbin
         write(14,'(f12.6,2x,g16.8)') (i-0.5_wp)*DL1,RDF(i)/(natom-nseed)
      end do
      write(io,'(/)')
      call flush(io)
   END SUBROUTINE RDF_PRINT

END MODULE RDF_MOD

