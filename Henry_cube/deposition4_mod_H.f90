
MODULE deposition_mod
   USE Henrys_law_calc_mod
   USE precision_mod, only: wp
   implicit none
   logical,parameter:: use_voidage_calc = .TRUE.
   integer,parameter,private:: ntrial = 10000
   integer,parameter,private:: nbinmax = 14
   integer,parameter:: iOSiOH3 = 6
   integer,parameter:: iOSiO = 7
   integer,parameter:: iOHB = 8
   integer:: nbin,itop_nonp,itop_nonp_SiO4,itop_nonp_H2O
   real(wp):: nden(nbinmax,0:4),voidage(nbinmax)
   real(wp):: ndenrxn(nbinmax,0:8)
   real(wp):: henryk(nbinmax,nbinmax,nbinmax,4),inthenryk(nbinmax,4),volbin
   real(wp):: intvoidage(nbinmax)
CONTAINS

   elemental function voidage_approx(no,sigmap)
      real(wp):: voidage_approx
      real(wp),intent(in):: no,sigmap
      real(wp),parameter:: alpha = 34.7_wp
      real(wp),parameter:: ro = 1.5_wp
      real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
      real(wp):: t1,t2
      t1 = -(4.0_wp/3.0_wp)*pi*no*(ro*ro*ro)
      t2 = 1.0_wp + alpha*no*(sigmap/(2.0_wp*ro))
      voidage_approx = exp(t1*t2*t2*t2)
   end function voidage_approx


   SUBROUTINE den_profile(ibegin,iend,dl)
      USE coordinates_mod
      USE atom_types
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: dl
      real(wp):: top_atom,delzi,delz
      integer:: j,ii,jj
      top_atom = maxval(rxyz(ibegin:iend,3))
      nbin = int(top_atom/dl) + 1
      delz = dl
      delzi = 1.0_wp/delz
      nden(1:nbin,:) = 0.0_wp
      do j = ibegin,iend
         ii = int(rxyz(j,3)*delzi + 1.0_wp)
         if (ii > nbin) then
            if ( rxyz(j,3) == top_atom) then
               ii = nbin
            else
               stop 'ii > nbin'
            end if
         end if
         jj = atom(j)
         nden(ii,jj) = nden(ii,jj) + 1.0_wp
      end do
      nden(1:nbin,:) = nden(1:nbin,:)/(boxl*boxl*delz*angstrom**3)
   END SUBROUTINE den_profile


   SUBROUTINE integrate_fn(ibegin,iend,fv,ifv)
      integer,intent(in):: ibegin,iend
      real(wp),intent(in):: fv(:)
      real(wp),intent(out):: ifv(:)
      real(wp):: tot
      integer:: i
      tot = 0.0_wp
      do i = ibegin,iend
         tot = tot + fv(i)
         ifv(i) = tot
      end do
      ifv(ibegin:iend)=ifv(ibegin:iend)/tot
   END SUBROUTINE


!   SUBROUTINE voidage_profile(probe,dl)


    SUBROUTINE henry_profile(probe,dl,ip)
      USE voidage_mod
      USE coordinates_mod
      USE probe_mol_mod
      USE nlist_mod
      type(probe_mol),intent(in):: probe
      real(wp),intent(in):: dl
      integer,intent(in):: ip
      real(wp):: facc,uljel,khenry,fover
      real(wp):: xl,xu,yl,yu,zl,zu
      integer:: i,j,k,ntrial,nbinx,nbiny,nbinz
      nbinx = ncelx
      nbiny = ncely
      nbinz = ncelz
      ntrial = 10000
      fover = 0.8_wp
      do k = 1,nbinz
         zl = (k-1)*delz
         zu = k*delz
         do j = 1,nbiny
            yl = (j-1)*dely-boxl2
            yu = j*dely-boxl2
            do i = 1,nbinx
               xl = (i-1)*delx-boxl2
               xu = i*delx-boxl2
               call Henrys_law_calc(ntrial,xl,xu,yl,yu,zl,zu,probe,fover,Uljel,Khenry,facc)
               henryk(i,j,k,ip) = Khenry
print '(3i5,3f14.8)',k,j,i,Khenry,Uljel,facc
            end do
         end do
      end do
   END SUBROUTINE



   SUBROUTINE prob_dist(fv,intfv,void_crit,itop_nonp)
      real(wp),intent(in):: fv(:)
      real(wp),intent(inout)::intfv(:)
      real(wp),intent(in):: void_crit
      integer,intent(out):: itop_nonp
      integer:: i
      itop_nonp = 1
      do i = nbin,1,-1
         if ( fv(i) < void_crit ) then
            itop_nonp = i+1
            exit
         end if
      end do
      call integrate_fn(itop_nonp,nbin,fv,intfv)
   END SUBROUTINE


   SUBROUTINE sample_z_bin(itop_nonp,iz)
      USE rand_mod
      integer,intent(in):: itop_nonp
      integer,intent(out):: iz
      integer:: iv1(1)
      real(wp):: rnd
      rnd = rand()
      iv1 = MINLOC(ABS(intvoidage(itop_nonp:nbin)-rnd))+itop_nonp-1
      iz = iv1(1)
   END SUBROUTINE

   SUBROUTINE sample_z_bin_uniform(itop_nonp,iz)
      USE rand_mod
      integer,intent(in):: itop_nonp
      integer,intent(out):: iz
      real(wp):: rnd
      rnd = rand()
      iz = int(rnd*(nbin - itop_nonp + 1)) + itop_nonp
   END SUBROUTINE

END MODULE deposition_mod

