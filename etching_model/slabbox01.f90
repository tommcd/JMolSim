
include 'precision_mod.f90'
include 'command_line_mod.f90'

!!>include 'command_arg_mod.f90'

MODULE command_arg_mod
   USE precision_mod
   USE command_line_mod
   implicit none
   integer, private:: icur_arg = -1
   !
   ! Gets a command argument and assigns it to a variable.
   ! It checks the status of the argument and writes a record to stdout.
   !
   public:: next_command_argument
   interface next_command_argument
      module procedure next_command_argument_int, next_command_argument_real, &
                       next_command_argument_logical,next_command_argument_string
   end interface
!
CONTAINS

   SUBROUTINE next_command_argument_int(var,var_name,stat)
      integer,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_int

   SUBROUTINE next_command_argument_real(var,var_name,stat)
      real(wp),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_real

   SUBROUTINE next_command_argument_logical(var,var_name,stat)
      logical,intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
         stop
      end if
      read(unit=carg,fmt=*) var
      write (*,*) var_name,' = ',var
   END SUBROUTINE next_command_argument_logical

   SUBROUTINE next_command_argument_string(var,var_name,stat)
      character(*),intent(out):: var
      integer,intent(out):: stat
      character(len=*),intent(in):: var_name
      character(len=132):: carg
      integer:: nc
      icur_arg = icur_arg + 1
      call get_command_argument(icur_arg, carg, length = nc, status = stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', icur_arg
         write (*,*) 'arg = ', carg(1:nc)
      end if
      var = ''
      var(1:nc) = carg(1:nc)
      write (*,*) var_name,' = ',trim(var)
   END SUBROUTINE next_command_argument_string

END MODULE command_arg_mod

include 'Keating_parameters_mod_vonAlfthan.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'rotate_axis_mod.f90'
include 'tetra_coords_mod.f90'
include 'global_vars_mod.f90'
include 'HKNonLattice2.f90'
include 'seaton_mod2.f90'
include 'atom_types_mod.f90'
include 'constants_mod2.f90'

!include 'coordinates_mod2.f90'

MODULE coordinates_mod
   USE precision_mod, only: wp
   USE global_vars_mod, only: natom,angstrom
   implicit none
   real(wp),allocatable:: rxyz(:,:),vxyz(:,:),fxyz(:,:)
   real(wp):: boxl,boxli,boxl2
   interface pbc
      module procedure  pbc_v, pbc_a, pbc_n_v, pbc_n_a
   end interface
CONTAINS

   PURE SUBROUTINE pbc_v(r)
      real(wp),intent(inout):: r(:)
      r(1:3) = r(1:3) - boxl*anint(r(1:3)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_a(r)
      real(wp),intent(inout):: r(:,:)
      r(:,1:3) = r(:,1:3) - boxl*anint(r(:,1:3)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_n_v(r,n)
      integer,intent(in):: n
      real(wp),intent(inout):: r(:)
      r(1:n) = r(1:n) - boxl*anint(r(1:n)*boxli)
   END SUBROUTINE

   PURE SUBROUTINE pbc_n_a(r,n)
      integer,intent(in):: n
      real(wp),intent(inout):: r(:,:)
      r(:,1:n) = r(:,1:n) - boxl*anint(r(:,1:n)*boxli)
   END SUBROUTINE

   SUBROUTINE WRITE_XYZ(iu,nattached)
      USE atom_types_mod, only: atom_name,atom
      integer,intent(in):: iu,nattached
      integer:: i
      write(iu,*) natom
      write(iu,'(a,i6)') 'Silica   nattached = ',nattached
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)),rxyz(i,1:3)*Angstrom
      end do
      write(iu,*)
      call flush(iu)
   END SUBROUTINE WRITE_XYZ

   SUBROUTINE READ_XYZ(iu,nattached)
      integer,intent(in):: iu,nattached
      integer:: i,it
      real(wp):: x,y,z
      character(32):: ctmp
      read(iu,*) natom
      read(iu,*) ctmp,ctmp,ctmp,it
      if (it /= nattached) stop 'it /= nattached'
      do i = 1,natom
         read(iu,*) ctmp,x,y,z
         rxyz(i,1:3) = (/ x,y,z /)/Angstrom
      end do
      rewind(iu)
   END SUBROUTINE READ_XYZ

END MODULE coordinates_mod


include 'connectivity_mod2.f90'
include 'charges_implicitH_mod.f90'
include 'probe_mol_mod2.f90'
include 'probe_mol_init_mod.f90'
include 'atom_list_mod.f90'
include 'bond_angle_types_mod.f90'
include 'bond_angle_list_mod.f90'
include 'frames_mod2.f90'

!!>include 'force_keating.f90'

MODULE force_keating_mod
   USE precision_mod
   USE constants_mod
   implicit none
CONTAINS

   SUBROUTINE energy_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: rl1(3),r1sq,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         rl1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         ebond = ebond + kbond(ib)*(r1sq - abond(ib))**2
      end do
      ebond = 0.5_wp*ebond
!
!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         eangle = eangle + kang(ia)*(cos_theta - ctheta(ia))**2
      end do
      eangle = 0.5_wp*eangle
   END SUBROUTINE energy_keating


   SUBROUTINE force_keating(ebond,eangle)
      USE bond_angle_list_mod
      USE Keating_parameters_mod
      USE coordinates_mod
      real(wp),intent(out):: ebond,eangle
      integer:: ib,ia,ia1,ia2,ia3
      real(wp):: k,dist,rl1(3),r1sq,ktheta,rl2(3),r2sq
      real(wp):: rl1l2,cos_theta,coef,coef1,coef2,coef3
!
!     Keating bond streching
      ebond = 0.0_wp
      do ib = 1,nbondtot
         k = KSiO
         dist = ASiO
         ia1 = ibond(1,ib)
         ia2 = ibond(2,ib)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         r1sq = sqrt(dot_product(rl1,rl1))
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - k*(r1sq - dist)*rl1(1:3)/r1sq
         fxyz(ia2,1:3) = fxyz(ia2,1:3) + k*(r1sq - dist)*rl1(1:3)/r1sq
         ebond = ebond + 0.5_wp*k*(r1sq - dist)**2
      end do

!     Keating bond angle bending
      eangle = 0.0_wp
      do ia = 1,nang
         ktheta = kang(ia)
         ia1 = iang(1,ia)
         ia2 = iang(2,ia)
         ia3 = iang(3,ia)
         rl1 = rxyz(ia1,1:3) - rxyz(ia2,1:3)
         rl2 = rxyz(ia3,1:3) - rxyz(ia2,1:3)
         call pbc(rl1)
         call pbc(rl2)
         r1sq = sqrt(dot_product(rl1,rl1))
         r2sq = sqrt(dot_product(rl2,rl2))
         rl1l2 = dot_product(rl1,rl2)
         cos_theta = (rl1l2)/(r1sq*r2sq)
         coef = ktheta*(cos_theta - ctheta(ia))
         coef1 = r1sq**3*r2sq
         coef2 = r1sq*r2sq
         coef3 = r1sq*r2sq**3
         fxyz(ia2,1:3) = fxyz(ia2,1:3) - coef*(-(rl1 + rl2)/coef2 &
                                            + rl1l2*rl1/coef1 &
                                            + rl1l2*rl2/coef3)
         fxyz(ia1,1:3) = fxyz(ia1,1:3) - coef*(rl2/coef2 - rl1l2*rl1/coef1)
         fxyz(ia3,1:3) = fxyz(ia3,1:3) - coef*(rl1/coef2 - rl1l2*rl2/coef3)
         eangle = eangle + 0.5_wp*ktheta*(cos_theta - ctheta(ia))**2
      end do
   END SUBROUTINE force_keating

END MODULE force_keating_mod


include 'nlist_mod2.f90'
include 'repul_force_NL_vonAlfthan2.f90'
!include 'quat2mat_mod.f90'
!include 'ran_point_sphere_mod.f90'
!include 'insert_mod.f90'
!include 'lj_el_mod_O2_N2_CO2.f90'
!include 'force_en_el_mod.f90'



MODULE condensation_mod
   USE precision_mod, only: wp
   implicit none
CONTAINS

   SUBROUTINE condensation(iO1,iO2)
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      integer,intent(inout):: iO1,iO2
      integer:: i,ii,ilst(3),itmp
      integer:: iSi2,iSi1
!
!      H          H        O
!      |    +     |  -->  / \  + H2O
!   Si-O       Si-O
!
!
!      |        |           iO1           iO2
!     iO1   +  iO2   -->    /  \   +      /  \
!     /        /           /    \
!   iSi1     iSi2        iSi1  iSi2
!
!! iO2 cannot be a seed
!      if (iO2 <= nseed) then
!         if (iO1 <= nseed) RETURN  ! don't allow 2 seed OH groups to react
!         itmp = iO2
!         iO2 = iO1
!         iO1 = itmp
!      end if

!     Pick iSi1, the Silicon attached to iO1
      iSi1 = 0
      do i = 1,ncmax(atom(iO1))
         ii = proximity(iO1,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi1 = ii
         end if
      end do
      iSi2 = 0
      do i = 1,ncmax(atom(iO2))
         ii = proximity(iO2,i)
         if (ii == 0) cycle
         if (atom(ii) == iSilicon) then
            iSi2 = ii
         end if
      end do
!
      atom(iO1) = iOxygen
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO1, 0, iSi2)
      call delete_atom(iO2)
      RETURN
   END SUBROUTINE

END MODULE condensation_mod


MODULE bond_switch_O_mod
   USE precision_mod, only: wp
   IMPLICIT NONE
CONTAINS

   SUBROUTINE bond_switch_O(iOb,success)
      USE atom_types_mod
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      implicit none
      integer,intent(in):: iOb
      logical,intent(out):: success
      integer:: iO1,iO2
      integer:: iSi1,iSi2,iSi3,iSi4
!     logical:: consistent

! TGV method
! Bond switch for SiO2 proposed by:
! Yuhai Tu, G. Grinstein, and D. Vanderbilt, Phys. Rev. Lett. 81,4899 (1998).
!
!      O        O2              O  O2
!      |        |               | /
!  O--Si1--Ob--Si2--O  -->  O--Si1--Ob--Si2--O
!      |        |                      / |
!      O1       O                    O1  O
!

      iSi1 = proximity(iOb,1)
      iSi2 = proximity(iOb,2)
      if (iSi1 == 0) stop 'oxygen malfunction 1'
      if (iSi2 == 0) stop 'oxygen malfunction 2'

! randomly select an oxygen iO1 attached to iSi1 (not iOb)
      do
         iO1 = proximity(iSi1,int(rand()*ncmax(atom(iSi1))) + 1)
         if (iO1 /= iOb) exit
      end do
! repeat for iSi2
      do
         iO2 = proximity(iSi2,int(rand()*ncmax(atom(iSi2))) + 1)
         if (iO2 /= iOb) exit
      end do

! no 2-Si rings are allowed to form
      iSi3 = proximity(iO2,1)
      if (iSi3 == iSi2) iSi3 = proximity(iO2,2)
      if (nearest_neighbor2(iSi3,iSi1)) then
         success = .FALSE.
         RETURN
      end if
      iSi4 = proximity(iO1,1)
      if (iSi4 == iSi1) iSi4 = proximity(iO1,2)
      if (nearest_neighbor2(iSi4,iSi2)) then
         success = .FALSE.
         RETURN
      end if

      call set_proximity(iSi1,iO1,iO2)
      call set_proximity(iO1,iSi1,iSi2)
      call set_proximity(iSi2,iO2,iO1)
      call set_proximity(iO2,iSi2,iSi1)
      success = .TRUE.

   END SUBROUTINE bond_switch_O

END MODULE bond_switch_O_mod



PROGRAM SIMBOX
      USE precision_mod
      USE command_arg_mod
      USE atom_types_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE nlist_mod
      USE bond_angle_list_mod
      USE rand_mod
      USE HKNonLattice_mod
      USE frames_mod
      USE condensation_mod
      USE bond_switch_O_mod
      USE probe_mol_mod
      USE probe_mol_init_mod
      implicit none
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: NAvo = 6.02214e23_wp
      real(wp):: n_si,n_ob,n_oh
      real(wp):: r3(3),rSiO(3),rn1,rn2,rn3
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      integer:: natom_old,natom_old1
      integer,allocatable:: atom_old(:),proximity_old(:,:),proximity2(:,:)
      integer,allocatable:: atom_old1(:),proximity_old1(:,:)
      real(wp),allocatable:: rxyz_old(:,:),bondl(:),rxyz_old1(:,:)
      integer,allocatable:: atomc(:)
      integer,allocatable:: cluster_list_x(:),cluster_list_y(:),cluster_list_z(:)
      real(wp):: CL,bl1,bl2,r3l,ztrans,dens_gcm3,volbox,fdens_wt,bondlmax
      real(wp):: r1(3),r2(3),rmax,rmax2,rl1,rl2,rmaxSiSi,rmaxSiSi2
      integer:: narg,stat,noxdel,iSi1,iSi2,iSi,iOx
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz,isil,jj,ibmax(1)
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z
      integer:: itmp,ibx,icx,iby,icy,ibz,icz,l,natcl,lst4(4),nc
      integer:: ipairs(4),ip(3,4),itry(3),ntry,ia,i1,i2,a1,a2
      integer:: icx_n,icy_n,icz_n,icl,icly,iclz,irand_seed0
      integer:: n_oh_top, n_oh_bot, n_Si_top(4), n_Si_bot(4), n_Si_OH_top(0:4), n_Si_OH_bot(0:4)
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster,success
      character(len=132) infile,ctmp,outfile,c6*6,c5*5,ftyp*3,carg
      type(atom_list):: lstSi
!
      ip(1,:) = (/ 1,2,3,4 /)
      ip(2,:) = (/ 1,4,3,2 /)
      ip(3,:) = (/ 1,3,2,4 /)
!
      narg = command_argument_count()
      write (*,*) 'number of command arguments = ', narg
      call next_command_argument(carg, "command",stat)
      if (narg /= 6) then
         write(*,*)'usage :'
         write(*,*) trim(carg),' ztrans rmax rmaxSiSi irand_seed  infile  outfile'
         write(*,*)
         write(*,*) "ztrans = translate z-dir"
         write(*,*) "rmax = maximum radius for condensing OH's [Angstrom]"
         write(*,*) 'max radius for Silicas [Angstrom]'
         write(*,*) 'irand_seed = seed for RNG'
         write(*,*) 'infile = name of input file'
         write(*,*) 'outfile = base name of output file'
         stop 'error: incorrect no. of command line arguments'
      end if
      call next_command_argument(ztrans, "translate z-dir",stat)
      call next_command_argument(rmax, "max radius for condensing OH's",stat)
      call next_command_argument(rmaxSiSi, "max radius for Silicas",stat)
      call next_command_argument(irand_seed0, "seed for RNG",stat)
      call next_command_argument(infile, "input file",stat)
      call next_command_argument(outfile, "base name of output file",stat)
      irand_seed = irand_seed0
      write(*,'(a,2f9.6)') 'Check: first 2 random numbers are ',rand(),rand()
      rmax = rmax/Angstrom
      rmax2 = rmax**2
      rmaxSiSi = rmaxSiSi/Angstrom
      rmaxSiSi2 = rmaxSiSi**2

!
      open(unit=14,file=trim(infile))
      nseed = 0
!
      nc = len_trim(infile)
print *,"'",trim(infile),"'"
print *,nc
      ftyp = infile(nc - 2:nc)
      select case(ftyp)
      case('xyz')
         read (14,*) natom
         read (14,*) boxl
         open(unit=15,file=trim(infile(1:nc - 3))//'con')
      case('pdb')
         read (14,*) c6,natom  ;print *,'natom = ',natom
         read (14,*) c6,boxl   ;print *,' boxl = ',boxl
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
!
      natom_max = natom + 1000
      boxl = boxl/Angstrom
      boxl2 = boxl*0.5_wp
      boxli = 1.0_wp/boxl
      volbox = (boxl*boxl*boxl)*(Angstrom**3)
!
      allocate(rxyz(1:natom_max,3),rxyz_old(1:natom_max,3))
      allocate(atom(1:natom_max),atom_old(1:natom_max))
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(cluster_list_x(natom_max),cluster_list_y(natom_max),cluster_list_z(natom_max))
      allocate(proximity(natom_max,4),proximity_old(natom_max,4),proximity2(natom_max,4))
      allocate(rxyz_old1(1:natom_max,3))
      allocate(atom_old1(1:natom_max))
      allocate(proximity_old1(natom_max,4))
!
      select case(ftyp)
      case('xyz')
         do i = 1,natom
            read (14,*) ctmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
         close(14)
      case('pdb')
         do i = 1,natom
            read (14,'(a6)',advance = 'no') c6
            read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
      end select
      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom

      proximity = 0
      select case(ftyp)
      case('xyz')
         read(15,*) itmp
         do i = 1,natom
            read(15,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(15)
      case('pdb')
         do i = 1,natom
            read(14,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(14)
      end select
!
! randomly translate the initial structure
!     rn1 = rand()
!     rn2 = rand()
!     rn3 = rand()
!     rxyz(1:natom,1) = rxyz(1:natom,1) + (2.0_wp*rn1 - 1.0_wp)*boxl
!     rxyz(1:natom,2) = rxyz(1:natom,2) + (2.0_wp*rn2 - 1.0_wp)*boxl
!     rxyz(1:natom,3) = rxyz(1:natom,3) + (2.0_wp*rn3 - 1.0_wp)*boxl
!     call pbc(rxyz)
!
! translate by ztrans in z-direction
      rxyz(1:natom,3) = rxyz(1:natom,3) + ztrans
      call pbc(rxyz)


      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
!
! count the clusters
      call Init_HKNonLattice(natom_max)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      print *,'n_cluster = ',n_cluster
      if (n_cluster /= 1) then
         stop 'initially there should be one cluster'
      end if

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
      call check_proximity2(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent (2) ',ii
         stop
      end if

!do i = 1,natom
!if (atom(i) ==iSilicon) then
!   do j =1,4
!      ii = proximity(i,j)
!      if (ii == 0) STOP '0 error in proximity'
!   end do
!end if
!end do
!do i = 1,natom
!if ( all(proximity(i,:) == 0) ) then
!   write(*,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
!   stop 'error in proximity'
!end if
!end do


      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom,natom_max)
      call NEW_NLIST(1,natom)

! Make array of silicon atoms
      lstSi%n = 0
      do i = 1,natom
         if (atom(i) == isilicon) call add_to_list(lstSi,i)
      end do
      print *,'no. of silicons = ',lstSi%n
      print *,count(atom(1:natom) == isilicon)
!!
!! Perfom the etching procedure
!!
!      isil = 0
!      etching: do
!         isil = isil + 1
!
!         ! Make array of silicon atoms
!         lstSi%n = 0
!         do i = 1,natom
!            if (atom(i) == isilicon) call add_to_list(lstSi,i)
!         end do
!
!         ! pick a Silicon to remove
!100      CONTINUE
!         ! print *,'lstSi%n ',lstSi%n
!         if (lstSi%n == 0) STOP 'NO suitable Silicon found for removal'
!         call rand_from_list(lstSi,iat)
!!do
!!   iat = int(rand()*natom) + 1
!!   if (atom(iat) == iSilicon) exit
!!end do
!         lst4 = proximity(iat,:)
!
!! store old config
!         call store_config()
!
!         where (lst4 == natom) lst4 = iat
!         call delete_atom(iat)
!
!         itry = (/ 1,2,3 /)
!         ntry = 3
!         success = .false.
!         do ii = 1,3
!            jj = int(rand()*ntry) + 1
!            ipairs = lst4(ip(itry(jj),:))
!            RING_if: if (OH_groups_OK(ipairs(1),ipairs(2)).and.OH_groups_OK(ipairs(3),ipairs(4))) then
!               ! Check if the 2 pairs of OH are are close enough
!               r1 = rxyz(ipairs(1),:) - rxyz(ipairs(2),:)
!               r2 = rxyz(ipairs(3),:) - rxyz(ipairs(4),:)
!               call pbc(r1)
!               call pbc(r2)
!               rl1 = dot_product(r1,r1)
!               rl2 = dot_product(r2,r2)
!               OH_if: if ( max(rl1,rl2) < rmax2 ) then
!                  ! Check if the 1st two Si are are close enough
!                  iSi1 = 0
!                  do i = 1,2
!                     k = proximity(ipairs(1),i)
!                     if (k == 0) cycle
!                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
!                     iSi1 = k
!                  end do
!                  iSi2 = 0
!                  do i = 1,2
!                     k = proximity(ipairs(2),i)
!                     if (k == 0) cycle
!                     ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
!                     iSi2 = k
!                  end do
!                  r1 = rxyz(iSi1,:) - rxyz(iSi2,:)
!                  call pbc(r1)
!                  rl1 = dot_product(r1,r1)
!                  Si_1_if: if ( rl1 < rmaxSiSi2 ) then
!                     ! Check if the 2nd two Si are are close enough
!                     iSi1 = 0
!                     do i = 1,2
!                        k = proximity(ipairs(3),i)
!                        if (k == 0) cycle
!                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
!                        iSi1 = k
!                     end do
!                     iSi2 = 0
!                     do i = 1,2
!                        k = proximity(ipairs(4),i)
!                        if (k == 0) cycle
!                        ! if (atom(k) /= iSilicon) STOP 'atom(k) == iSilicon'
!                        iSi2 = k
!                     end do
!                     r2 = rxyz(iSi1,:) - rxyz(iSi2,:)
!                     call pbc(r2)
!                     rl2 = dot_product(r2,r2)
!                     Si_2_if: if ( rl2 < rmaxSiSi2 ) then
!                        success = .true.
!                        exit
!                     else Si_2_if
!                        ! print *,'Si - Si 2 failed'
!                     end if Si_2_if
!                  else Si_1_if
!                     ! print *,'Si - Si 1 failed'
!                  end if Si_1_if
!               else OH_if
!                  ! print *,'OH constraint failed'
!               end if OH_if
!            else RING_if
!               ! print *,'Ring constraint failed'
!               itry(jj) = itry(ntry)
!               ntry = ntry - 1
!            end if RING_if
!         end do
!         if (.not.success) then
!            call restore_config()
!            call remove_from_list(lstSi,iat)
!            GOTO 100
!!            print *,'could not find 2 suitable pairs to condense'
!!            cycle etching
!         end if
!
!         r1 = rxyz(ipairs(2),:) - rxyz(ipairs(1),:)
!         call pbc(r1)
!         rxyz(ipairs(1),:) = rxyz(ipairs(1),:) +  r1*0.5_wp
!         call pbc(rxyz(ipairs(1),:))
!
!         r2 = rxyz(ipairs(4),:) - rxyz(ipairs(3),:)
!         call pbc(r2)
!         rxyz(ipairs(3),:) = rxyz(ipairs(3),:) +  r2*0.5_wp
!         call pbc(rxyz(ipairs(3),:))
!
!         if (maxval(ipairs(1:2)) > maxval(ipairs(3:4))) then
!            call condensation(ipairs(1),ipairs(2))
!            call condensation(ipairs(3),ipairs(4))
!         else
!            call condensation(ipairs(3),ipairs(4))
!            call condensation(ipairs(1),ipairs(2))
!         end if
!
!
!         call bond_list
!         if (.not.allocated(bondl)) allocate(bondl(nbondtot))
!
!         do ib = 1,nbondtot
!            r1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
!            call pbc(r1)
!            bondl(ib) = sqrt(dot_product(r1,r1))
!            ! print *,ib,ibond(1,ib),ibond(2,ib), bondl(ib)
!            ! energy = 0.5_wp*KSiO*(r1sq - ASiO)**2
!         end do
!         ibmax = maxloc(bondl(1:nbondtot))
!!         print *,bondl(ibmax(1))
!         bondlmax = bondl(ibmax(1))
!
!         call store_config1()
!         call bond_switch_O(ipairs(1),success)
!         call bond_switch_O(ipairs(3),success)
!
!         call bond_angle_list()
!         call bond_list
!         if (.not.allocated(bondl)) allocate(bondl(nbondtot))
!
!         do ib = 1,nbondtot
!            r1 = rxyz(ibond(1,ib),1:3) - rxyz(ibond(2,ib),1:3)
!            call pbc(r1)
!            bondl(ib) = sqrt(dot_product(r1,r1))
!         end do
!         ibmax = maxloc(bondl(1:nbondtot))
!         call pbc(r1)
!         if (bondl(ibmax(1)) > bondlmax) then
!            call restore_config1()
!         else
!            print *,'switch accepted'
!         end if
!
!
!
!
!
!
!!call check_proximity(ok,ii)
!!if (.not.ok) then
!!print *,'just after condensation 2'
!!print *,'proximity array is inconsistent ',ii
!!stop
!!end if
!!do i = 1,natom
!!ia = atom(i)
!!select case(ia)
!!case(iSilicon)  ! Check Si is bonded to 4 atoms
!!   do j = 1,ncmax(ia)
!!      k = proximity(i,j)
!!      if (k == 0) STOP 'ERROR 01'
!!      if (atom(k) ==iSilicon) STOP 'ERROR 02'
!!      if (atom(k) ==iHydrogen) STOP 'ERROR 03'
!!   end do
!!case(iOxygen)  ! Bridging Oxygen must be bonded to Si
!!   i1 = proximity(i,1)
!!   if (i1 == 0) then
!!      print *,'i = ',i
!!      print *,'ipairs = ',ipairs
!!      STOP 'ERROR 04 B'
!!   end if
!!   a1 = atom(proximity(i,1))
!!   if (a1 /= iSilicon) STOP 'ERROR 05'
!!      i2 = proximity(i,2)  !  must be bonded to  2 Silicons
!!      if (i2 == 0) STOP 'ERROR 06'
!!      a2 = atom(proximity(i,2))
!!      if (a2 /= iSilicon) STOP 'ERROR 07'
!!   if (any(proximity(i,3:) /= 0)) STOP 'ERROR 08'
!!case default
!!   stop 'unknown atom type'
!!end select
!!end do
!
!
!!! remove the free oxygens
!!         j = natom
!!         noxdel = 0
!!         do
!!            if (atom(j) == iOxygen .or. atom(j) == iOxygenH) then
!!            if (count(proximity(j,:) > 0) < 1) then
!!               call delete_atom(j)
!!               noxdel = noxdel + 1
!!            end if
!!            end if
!!            j = j - 1
!!            if (j == 0) exit
!!         end do
!!! convert singly bonded oxygens to OH
!!         do i = 1,natom
!!            if (atom(i) == iOxygen) then
!!            if (count(proximity(i,:) > 0) < 2) then
!!               atom(i) = iOxygenH
!!            end if
!!            end if
!!         end do
!
!! check the clusters
!         atomL = 0
!         call Init_HKNonLattice(natom)
!         call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!         !print *,'n_cluster = ',n_cluster
!         if (n_cluster > 1) then
!!           call remove_from_list(lstSi,iat)
!!           if (atom(iat) ==iSilicon) then
!!              where(lstSi%i == (natom+1)) lstSi%i = iat
!!           end if
!            call restore_config()
!            cycle etching
!         end if
!
!         ! check for errors
!         if (ANY(proximity(1:natom,:) > natom)) then
!            print *,'error: proximity(1:natom,:) > natom '
!            do i = 1,natom
!               if (ANY(proximity(i,:) > natom)) then
!               print *,i,' : ',proximity(i,:)
!               end if
!            end do
!            stop
!         end if
!         if (ANY(proximity(natom + 1:,:) /= 0)) then
!            print *,'error: proximity(natom + 1:,:) /= 0 '
!            stop
!         end if
!
!
!      end do etching

!
! Use the Hoschen-Kopelman algorithm to label the clusters
!
      atomL = 0
      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
      print *,'n_cluster = ',n_cluster

      call bond_list
      proximity2 = proximity
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do i = 1,nbondtot
         r3 = rxyz(ibond(1,i),:) - rxyz(ibond(2,i),:)
         if (dot_product(r3,r3) > boxl2*boxl2 ) then
            do j = 1,4
               if (proximity2(ibond(1,i),j) == ibond(2,i)) then
                   proximity2(ibond(1,i),j) = 0
               end if
            end do
            do j = 1,4
               if (proximity2(ibond(2,i),j) == ibond(1,i)) then
                   proximity2(ibond(2,i),j) = 0
               end if
            end do
         end if
         if (ABS(r3(1)) > boxl2) then
            n_crossbond_x = n_crossbond_x + 1
            crossbond_x(n_crossbond_x) = i
         end if
         if (ABS(r3(2)) > boxl2) then
            n_crossbond_y = n_crossbond_y + 1
            crossbond_y(n_crossbond_y) = i
         end if
         if (ABS(r3(3)) > boxl2) then
            n_crossbond_z = n_crossbond_z + 1
            crossbond_z(n_crossbond_z) = i
         end if
      end do
      print *,'no. of crossbonds x,y,z: '
      print *,n_crossbond_x, n_crossbond_y, n_crossbond_z


      proximity2 = proximity
      n_crossbond_x = 0
      n_crossbond_y = 0
      n_crossbond_z = 0
      do ib = 1,nbondtot
         ii = ibond(1,ib)
         jj = ibond(2,ib)
         r3 = rxyz(ii,:) - rxyz(jj,:)
         call pbc(r3,2)
         if (dot_product(r3,r3) > boxl2*boxl2 ) then
            if (atom(ii)==iSilicon)then
               iSi = ii; iOx = jj
            else
               iSi = jj; iOx = ii
            end if
            where( proximity2(ii,:)==jj ) proximity2(ii,:) = 0
            where( proximity2(jj,:)==ii ) proximity2(jj,:) = 0
            natom = natom + 1
            if (natom > natom_max) stop 'natom_max too small'
            atom(natom) = atom(iOx)
            call change_connectivity(natom,0,iSi,proximity2)
            call change_connectivity(iSi,0,natom,proximity2)
            rSiO = rxyz(iOx,:) - rxyz(iSi,:)
            call pbc(rSiO)
            rxyz(natom,:) = rxyz(iSi,:) + rSiO
!
!            do k = 1,4
!               if (proximity2(ii,k) == jj) then
!                   proximity2(ii,k) = 0
!               end if
!            end do
!            do k = 1,4
!               if (proximity2(jj,k) == ii) then
!                   proximity2(jj,k) = 0
!               end if
!            end do
         end if
         if (ABS(r3(1)) > boxl2) then
            n_crossbond_x = n_crossbond_x + 1
            crossbond_x(n_crossbond_x) = ib
         end if
         if (ABS(r3(2)) > boxl2) then
            n_crossbond_y = n_crossbond_y + 1
            crossbond_y(n_crossbond_y) = ib
         end if
         if (ABS(r3(3)) > boxl2) then
            n_crossbond_z = n_crossbond_z + 1
            crossbond_z(n_crossbond_z) = ib
         end if
      end do
      print *,'no. of crossbonds x,y,z: '
      print *,n_crossbond_x, n_crossbond_y, n_crossbond_z

      proximity = proximity2



      atomL = 0

      call Init_HKNonLattice(natom)
      call HKNonLattice(natom,proximity2(1:natom,:),n_cluster,atomL)

      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
spanning_cluster = .true.
!!
!! Finding a cluster which crosses all faces
!!
!       icx_n = 0
!       icy_n = 0
!       icz_n = 0
!       do ibx = 1,n_crossbond_x
!          icx = crossbond_x(ibx)
!          if ( atomL(ibond(1,icx)) == atomL(ibond(2,icx)) ) then
!             icx_n = icx_n + 1
!             cluster_list_x(icx_n) = atomL(ibond(1,icx))
!          end if
!       end do
!
!       do iby = 1,n_crossbond_y
!          icy = crossbond_y(iby)
!          if ( atomL(ibond(1,icy)) == atomL(ibond(2,icy)) ) then
!             icy_n = icy_n + 1
!             cluster_list_y(icy_n) = atomL(ibond(1,icy))
!         end if
!       end do
!
!       do ibz = 1,n_crossbond_z
!          icz = crossbond_z(ibz)
!          if ( atomL(ibond(1,icz)) == atomL(ibond(2,icz)) ) then
!             icz_n = icz_n + 1
!             cluster_list_z(icz_n) = atomL(ibond(1,icz))
!          end if
!       end do
!
!      spanning_cluster = .false.
!      outer: do icl = 1, icx_n
!      do icly = 1, icy_n
!      do iclz = 1, icz_n
!         if ((cluster_list_x(icl) == cluster_list_y(icly)) .AND. &
!             (cluster_list_x(icl) == cluster_list_z(iclz))) then
!             spanning_cluster = .true.
!             exit outer
!         end if
!      end do
!      end do
!      end do outer
!
!      do icly = 1,natom
!         if ( atomL(icly) == cluster_list_x(icl) ) exit
!      end do
!
!      atomL = 0
!      call Init_HKNonLattice(natom)
!      call HKNonLattice(natom,proximity(1:natom,:),n_cluster,atomL)
!      print *,'n_cluster = ',n_cluster
!
!      j = natom
!      do
!         if (atomL(j) /= atomL(icly)) then
!print *,'deleting atoms not in 1st spanning cluster'
!            call delete_atom(j)
!print *,'delete ',j
!            atomL(j) = atomL(icly)
!         end if
!         j = j - 1
!         if (j == 0) exit
!      end do
!
!      do i = 1,natom
!         if ((atom(i) == iSilicon).and.(atomL(i) == atomL(icly))) then
!         do j = 1,4
!            ii = proximity(i,j)
!            if (ii == 0) STOP '5 error in proximity'
!         end do
!         end if
!      end do


! convert singly bonded oxygens to oh
      n_oh_top = 0
      n_oh_bot = 0
      n_Si_top = 0
      n_Si_bot = 0
      do i = 1,natom
         if (atom(i) == ioxygen) then
            if (count(proximity(i,:) > 0) < 2) then
               atom(i) = iOxygenH
               if(rxyz(i,3) > 0.0_wp)then
                  n_oh_top = n_oh_top + 1
               else
                  n_oh_bot = n_oh_bot + 1
               end if
            end if
         else if (atom(i) == iSilicon) then
            nc = count( proximity(i,:)==0 )
            if (count(proximity(i,:) > 0) < 4) then
               !atom(i) = 
               if(rxyz(i,3) > 0.0_wp)then
                  n_Si_top(nc) = n_Si_top(nc) + 1
               else
                  n_Si_bot(nc) = n_Si_bot(nc) + 1
               end if
            end if
         end if
      end do
      write(*,'(a,4i6)') 'n_Si_top ',n_Si_top
      write(*,'(a,4i6)') 'n_Si_bot ',n_Si_bot
      write(*,*) 'n_oh_top ',n_oh_top
      write(*,*) 'n_oh_bot ',n_oh_bot
      write(*,'(a,f12.6,a)') 'n_oh_top =',n_oh_top/(boxl*boxl),' nm^-2'
      write(*,'(a,f12.6,a)') 'n_oh_bot =',n_oh_bot/(boxl*boxl),' nm^-2'

      n_Si_OH_top = 0
      n_Si_OH_bot = 0
      do i = 1,natom
         if (atom(i) == iSilicon) then
            nc = noh(i)
            if(rxyz(i,3) > 0.0_wp)then
               n_Si_OH_top(nc) = n_Si_OH_top(nc) + 1
            else
               n_Si_OH_bot(nc) = n_Si_OH_bot(nc) + 1
            end if
         end if
      end do
      write(*,'(a,5i6)') '       # OH:',0,1,2,3,4
      write(*,'(a,5i6)') 'n_Si_OH_top ',n_Si_OH_top
      write(*,'(a,5i6)') 'n_Si_OH_bot ',n_Si_OH_bot

      write(*,*) 'natom ',natom
      write(*,*) 'irand_seed ',irand_seed
      if (spanning_cluster) then
         imve = 77
         write(unit=ctmp,fmt='(i6.6)') irand_seed0
         write(unit=c6,fmt='(i6.6)') natom
         write(unit=c5,fmt='(f5.2)') ztrans
         open(unit=imve,file=trim(outfile)//'_'//trim(adjustl(c5))//'_'//trim(c6)//'_'//trim(ctmp)//'.pdb')
         call write_frame(imve,1,natom)
      else
         print *,'A spanning cluster was NOT generated'
      end if


! check everything is ok
! check each Si has 4 bonds
!      do i = 1,natom
!         if ((atom(i) == iSilicon)) then
!         do j = 1,4
!            ii = proximity(i,j)
!            if (ii == 0) STOP '5 error in proximity'
!         end do
!         end if
!      end do
      call check_proximity(ok,ii)
      if (.not.ok) then
         print *,'proximity array is inconsistent ',ii
         stop
      end if
!      call check_proximity2(ok,ii)
!      if (.not.ok) then
!         call write_proximity(6,1)
!         print *,'proximity array is inconsistent (2) ',ii
!         stop
!      end if
!
! calculate the density
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      n_si = n_si/volbox
      n_ob = n_ob/volbox
      n_oh = n_oh/volbox
      write(*,*) 'N_Si = ',n_si,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*1e24_wp/NAvo
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'
      STOP 'normal termination'

CONTAINS

   SUBROUTINE store_config()
      atom_old = atom
      rxyz_old = rxyz
      proximity_old = proximity
      natom_old = natom
   END SUBROUTINE

   SUBROUTINE restore_config()
      atom = atom_old
      rxyz = rxyz_old
      proximity = proximity_old
      natom = natom_old
   END SUBROUTINE

   SUBROUTINE store_config1()
      atom_old1 = atom
      rxyz_old1 = rxyz
      proximity_old1 = proximity
      natom_old1 = natom
   END SUBROUTINE

   SUBROUTINE restore_config1()
      atom = atom_old1
      rxyz = rxyz_old1
      proximity = proximity_old1
      natom = natom_old1
   END SUBROUTINE

END PROGRAM SIMBOX

