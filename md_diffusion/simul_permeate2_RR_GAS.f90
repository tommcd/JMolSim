include 'precision_mod.f90'
include 'commandline_mod.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'HKNonLattice2.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'Keating_parameters_vonAlfthan.f90'
include 'bond_list_mod.f90'
include 'box_frames_mod.f90'
include 'comvel_mod.f90'
include 'force_keating.f90'
include 'vverlet_mod.f90'
include 'charges2_mod.f90'
include 'nlist_mod.f90'
include 'repul_force_NL_vonAlfthan.f90'
include 'quat2mat_mod.f90'
include 'ran_point_sphere_mod.f90'
include 'insert_mod.f90'
include 'lj_el_mod_O2_N2_CO2.f90'
include 'force_en_el_mod.f90'
include 'msd_mod3.f90'
include 'vcf_mod3.f90'

PROGRAM SIMBOX
      USE commandline_mod
      USE precision_mod
      USE seaton_mod
      USE constants_mod
      USE global_vars
      USE coordinates_mod
      USE connectivity_mod
      USE atom_list_mod
      USE bond_list_mod
      USE rand_mod
      USE atom_types
      USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE keating_parameters
      USE nlist_mod
      USE repul_energy_mod
      USE charges_mod
      USE force_en_el_mod
      USE insert_mod
      USE msd_mod
      USE vcf_mod
      implicit none
      integer,parameter:: number_moves = 260000
      integer,parameter:: n_short = 4
      integer,parameter:: nequib = 10000
      integer,parameter:: nsamp = 10
      integer,parameter:: nprint_all = 1000
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,qj,sumqj,T_kelvin,T_kinetic,sumx,sumy,sumz
      real(wp),allocatable:: rxyzc(:,:),mol_vec(:,:)
      integer,allocatable:: atomc(:)
      real(wp):: dtreal,CL,time,EK,EK_gas,rr(3),ULJEL,dt_short,mtmp
      real(wp),allocatable:: massi(:,:),GAS_massi(:,:,:),mass(:,:),GAS_mass(:,:,:),rad_GAS(:)
      integer,allocatable:: mol_atoms(:)
      integer:: j,i,k,m,ii,ic = 0,natomc,istep,len,stat,nat
      integer:: itmp,i_short,narg,Ndf
      real(wp),allocatable:: coprxyz(:,:)
      real:: tp(0:10),sumtp(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5
      character(len=132):: infile,carg
      logical:: restart
!
      ntau_msd = 10000
      ntau_vcf = 10000
!
      call cpu_time(tp(0))
      open(unit=24,file='final.pdb')
      open(unit=17,file='energy.out')
      open(unit=18,file='temperature.out')
      open(unit=1003,file='msd_GAS.out')
      open(unit=1006,file='vcf_GAS.out')
!
      narg = command_argument_count()
      call get_command_argument (0, carg, len, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat
         stop
      end if
      write (*,*) 'command name = ', carg(1:len)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) carg(1:len),'  natom  atom_file dtreal[fs] T[K]  {xyz_restart  vel_restart} '
         stop
      end if
!
      call get_command_argument (1, carg, len, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', 1
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) natom_max
      print *,'natom = ',natom_max
!
      call get_command_argument (2, infile, len, stat)
      if (stat /= 0) then
         write (*,*) 'get_command_argument failed: status = ', stat, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, carg, len, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', 3
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) dtreal

!
      call get_command_argument (4, carg, len, stat)
      if (stat /= 0) then
         write (*,*) 'Getting command name failed with status = ', stat, ' arg = ', 4
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) T_kelvin
!
      call get_command_argument (5, infile, len, stat)
      if (stat /= 0) then
         write (*,*) 'get_command_argument failed: status = ', stat, ' arg = ', 5
         restart = .false.
      else
         restart = .true.
         write (*,*) 'restarting from ', infile(1:len)
         open(unit=34,file=trim(infile(1:len)))
         call get_command_argument (6, infile, len, stat)
         if (stat /= 0) then
            print *,'restart velocity file missing'
            stop
         else
            open(unit=35,file=trim(infile(1:len)))
         end if
      end if

!
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!rcutel = boxl2
!rcutel2 = boxl2**2
!
      n_GAS = 4
      n_atom_gas = 5
      allocate(mol_atoms(n_atom_gas))
      allocate(GAS_atom(1:n_GAS,n_atom_gas))
      allocate(GAS_gas_charge(n_GAS,n_atom_gas))
      allocate(GAS_mass(1:n_GAS,n_atom_gas,3),GAS_massi(1:n_GAS,n_atom_gas,3))
      allocate(rad_GAS(n_atom_gas))
      mol_atoms = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
!      mol_atoms = (/ iSilicon,iOxygenH,iOxygenH /)
      allocate( mol_vec(3,n_atom_GAS) )
!mol_vec(:,1) =0.0_wp
!mol_vec(:,2) = (/0.0_wp,0.0_wp,ASiO/)
!mol_vec(:,3) = (/ 0.0_wp, -(2.0_wp*((1.0_wp/3.0_wp)*ASiO))*sqrt(2.0_wp), (1.0_wp/3.0_wp)*ASiO /)
      call tetra_coords5(ASiO, mol_vec)
print '(3f12.6)',mol_vec
print *,mol_atoms
print '(a)',atom_name(mol_atoms)
!     Initialize the gas molecules
      do i = 1,n_GAS
         !GAS_gas_charge(i,:) =
         GAS_atom(i,:) = mol_atoms(1:n_atom_gas)
         do j = 1, n_atom_gas
            GAS_mass(i,j,:) = amass(mol_atoms(j))
         end do
      end do
      GAS_massi = 1.0_wp/GAS_mass

!print *,'GAS_mass',GAS_mass
!print *,'GAS_massi',GAS_massi

      rad_GAS = sigma_2(mol_atoms)

      natom = natom_max
!
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
      nat = 2
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(GAS_xyz(1:n_GAS,1:n_atom_gas,1:3))
      allocate(r_GAS_uc(n_GAS,n_atom_gas,3))
      allocate(dr_GAS(n_GAS,n_atom_gas,3))
      allocate(rcm_GAS_uc(n_GAS,3))
      allocate(vcm_GAS(n_GAS,3))
      allocate(Imige_GAS_xyz(1:n_GAS,1:n_atom_gas,1:3))
      allocate(fxyz(1:natom_max,3),vxyz(1:natom_max,3))
      allocate(GAS_fxyz(1:n_GAS,1:n_atom_gas,1:3))
      allocate(GAS_vxyz(1:n_GAS,1:n_atom_gas,1:3))
      allocate(atom(1:natom_max))
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
      allocate(mass(natom_max,3),massi(natom_max,3))

      proximity = 0

      do i = 1,natom
         read (14,'(a6)',advance = 'no') c6
         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
         if (atom(i) == iSilicon) mass(i,:) = mSi
         if (atom(i) == iOxygen)  mass(i,:) = mOx
         if (atom(i) == iOxygenH)  mass(i,:) = mOH
      end do
      massi = 1.0_wp/mass

      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
      do i = 1,natom
         read(14,'(a32)') ctmp
         do j = 1,ncmax(atom(i))
            c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
            read( unit=c5,fmt=* ) proximity(i,j)
         end do
      end do
      close(14)

      call assign_charge(1,natom)

      write(*,*) 'sum charge = ',sum(charge(1:natom))
      write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
      write(*,*) 'n Si = ',count(atom == iSilicon)
      write(*,*) 'n O  = ',count(atom == iOxygen)
      write(*,*) 'n OH = ',count(atom == iOxygenH)
      write(*,*) 'ntot = ',count(atom == iSilicon) + count(atom == iOxygen) + count(atom == iOxygenH)
!
      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom)
      call NEW_NLIST

! Insertion
      call insert_mol_GAS(n_GAS,n_atom_gas,mol_vec,GAS_xyz,100000,rad_GAS)
!print *
!do i = 1,n_gas
!do j = 1,n_atom_gas
!print '(3f12.6)',(GAS_xyz(i,j,k),k=1,3)
!end do
!end do
!stop
!
      if (restart) then
         call read_box_frame(34,nmol = n_GAS,print_all = .true.)
         call read_vel(35,nmol = n_GAS,print_all = .true.)
      end if
!
! Set up the uncorrected coords
      do i = 1, n_GAS
         r_GAS_uc(i,1,:) = GAS_xyz(i,1,:)
         do j = 2,n_atom_GAS
            rr = GAS_xyz(i,j,:) - GAS_xyz(i,1,:)
            call pbc(rr)
            r_GAS_uc(i,j,:) =  r_GAS_uc(i,1,:) + rr
         end do
      end do

      if (make_movie) then ! write out first frame of movie
         iframe = 0
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_GAS,print_all = .true.)
      end if
!========================================================================
!
      call LJ_INIT

      call comvel(natom,n_GAS,n_atom_GAS,T_kelvin,vxyz,GAS_vxyz,mass,GAS_mass)
print *,'net gas momentum',sum(GAS_vxyz(:,:,1)*GAS_mass(:,:,1)),&
                           sum(GAS_vxyz(:,:,2)*GAS_mass(:,:,2)),&
                           sum(GAS_vxyz(:,:,3)*GAS_mass(:,:,3))
      SUMX = 0.0_wp
      SUMY = 0.0_wp
      SUMZ = 0.0_wp
      do i = 1,n_GAS
      do j = 1,n_atom_GAS
         SUMX = SUMX + GAS_vxyz(i,j,1)*GAS_mass(i,j,1)
         SUMY = SUMY + GAS_vxyz(i,j,2)*GAS_mass(i,j,2)
         SUMZ = SUMZ + GAS_vxyz(i,j,3)*GAS_mass(i,j,3)
      end do
      end do

print *,'new sum x y z'
print *,SUMX
print *,SUMY
print *,SUMZ
!print *
!print *
!print '(2f)',1.0_wp/(GAS_massi)
!print *
!print '(2f)',GAS_vxyz
!print *
         EK = 0.0_wp
         do m = 1, natom
            EK = EK + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         EK_gas = 0.0_wp
         do m = 1, n_GAS
            do i = 1,n_atom_gas
               EK_gas = EK_gas + (0.5_wp/GAS_massi(m,i,1))*dot_product(GAS_vxyz(m,i,:),GAS_vxyz(m,i,:))
            end do
         end do
         Ndf = 3*(natom + n_atom_gas*n_GAS) - 3
         T_kinetic = (EK + EK_gas)*2.0_wp/(K_ev*Ndf)
print *,'T_kinetic = ',T_kinetic
         vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
         EK = 0.0_wp
         do m = 1, natom
            EK = EK + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         EK_gas = 0.0_wp
         do m = 1, n_GAS
            do i = 1,n_atom_gas
               EK_gas = EK_gas + (0.5_wp/GAS_massi(m,i,1))*dot_product(GAS_vxyz(m,i,:),GAS_vxyz(m,i,:))
            end do
         end do
         Ndf = 3*(natom + n_atom_gas*n_GAS) - 3
         T_kinetic = (EK + EK_gas)*2.0_wp/(K_ev*Ndf)
print *,'T_kinetic = ',T_kinetic

      !GAS_vxyz(1,1,:) = 0.01_wp !####################################################################TEST

      if (make_movie) then ! write out first frame of velocities
         iframe = 0
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_GAS,print_all = .false.)
      end if

      call set_bond_angle_lists()
      call msd_init(n_GAS)
      call vcf_init(n_GAS)

      time = 0.0_wp
      dt = dtreal*dt_fs
      dt_short = dt/n_short     !0.25_wp*dt_fs

      fxyz = 0.0_wp
      GAS_fxyz = 0.0_wp
      energy = 0.0_wp
!
!      call FORCE_ENERGY_LJ_EL_GAS(1,natom,1,n_GAS,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_GAS,3,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
! Electrostatic
!      call FORCE_ENERGY_EL(1,natom,1,n_GAS,3,GAS_gas_charge,GAS_xyz,GAS_fxyz,ULJEL)

!======================================================== main loop
      main_loop: do istep = 1,number_moves

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi

      GAS_vxyz = GAS_vxyz + 0.5_wp*dt*GAS_fxyz*GAS_massi
!print *,'main0: dt*GAS_fxyz*GAS_massi',dt*GAS_fxyz*GAS_massi


      fxyz = 0.0_wp
      GAS_fxyz = 0.0_wp
      call cpu_time(tp(1))
      call force_keating
      call repul_energy2(1,natom)
      call GAS_stretching(n_GAS,n_atom_GAS,GAS_xyz,ASiO,kSiO)
      call GAS_angle_bending(n_GAS,n_atom_GAS,GAS_xyz,ctheta = -1.0_wp/3.0_wp,ktheta = KOSiO)

!=============================== short loop ===
      short_loop: do i_short = 1, n_short
         coprxyz = rxyz
         call vverlet_a(dt_short,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)
!
         call gas_vverlet_a(dt_short,GAS_massi,GAS_xyz,GAS_vxyz,GAS_fxyz,dr_GAS)
         call pbc(GAS_xyz)
         r_GAS_uc = r_GAS_uc + dr_GAS

         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3))) then
               !! print *,'NEW LIST'
               call NEW_NLIST
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         GAS_fxyz = 0.0_wp
         energy = 0.0_wp
!call cpu_time(tp(4))
         call force_keating
!call cpu_time(tp(5)); sumtp(4) = sumtp(4) + tp(5)-tp(4)
         call repul_energy2(1,natom)
!call cpu_time(tp(6)); sumtp(5) = sumtp(5) + tp(6)-tp(5)
         call GAS_stretching(n_GAS,n_atom_GAS,GAS_xyz,ASiO,kSiO)
!call cpu_time(tp(7)); sumtp(6) = sumtp(6) + tp(7)-tp(6)
         call GAS_angle_bending(n_GAS,n_atom_GAS,GAS_xyz,ctheta = -1.0_wp/3.0_wp,ktheta = KOSiO)
!call cpu_time(tp(8)); sumtp(7) = sumtp(7) + tp(8)-tp(7)
!
         call vverlet_b(dt_short,massi,vxyz,fxyz)
         call gas_vverlet_b(dt_short,GAS_massi,GAS_vxyz,GAS_fxyz)
!call cpu_time(tp(9)); sumtp(8) = sumtp(8) + tp(9)-tp(8)
!============================end short loop ===
      end do short_loop
      call cpu_time(tp(2)); sumtp(1) = sumtp(1) + tp(2) - tp(1)

      fxyz = 0.0_wp
      GAS_fxyz = 0.0_wp
!
!     call FORCE_ENERGY_LJ_EL_GAS(1,natom,1,n_GAS,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_GAS,3,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
! Electrostatic
!      call FORCE_ENERGY_EL(1,natom,1,n_GAS,3,GAS_gas_charge,GAS_xyz,GAS_fxyz,ULJEL)
!
      call cpu_time(tp(3)); sumtp(2) = sumtp(2) + tp(3) - tp(2)

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi
      GAS_vxyz = GAS_vxyz + 0.5_wp*dt*GAS_fxyz*GAS_massi
!print *,'main1: dt*GAS_fxyz*GAS_massi',dt*GAS_fxyz*GAS_massi
      time = time + dtreal
      if (istep <= nequib) then
         ! print *,sum(fxyz(1:natom,:),dim=1)
         EK = 0.0_wp
         do m = 1, natom
            EK = EK + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         EK_gas = 0.0_wp
         do m = 1, n_GAS
            do i = 1,n_atom_gas
               EK_gas = EK_gas + (0.5_wp/GAS_massi(m,i,1))*dot_product(GAS_vxyz(m,i,:),GAS_vxyz(m,i,:))
            end do
         end do
         Ndf = 3*(natom + n_atom_gas*n_GAS) - 3
         T_kinetic = (EK + EK_gas)*2.0_wp/(K_ev*Ndf)

!print *,'GAS_vxyz 1',GAS_vxyz(1,1,:)
!print *,'GAS_vxyz 1',GAS_vxyz(1,2,:)
!print *,'T_kinetic = ',T_kinetic
         if (mod(istep,nsamp) == 0) then
            write(*,*) istep,tp(2) - tp(1),tp(3) - tp(2),tp(3) - tp(0)
            write(17,'(8g18.8)') time, (EK + EK_gas), energy + repulenergy*0.5_wp, (EK + EK_gas) + energy + repulenergy*0.5_wp
            write(18,'(8g18.8)') time, T_kinetic
            call flush(17); call flush(18)
         end if
!if (T_kinetic==0.0_wp) T_kinetic=1.0_wp
         vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
      end if

      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then
         do i = 1,n_GAS
            mtmp = 0.0_wp
            rcm_GAS_uc(i,:) = 0.0_wp
            vcm_GAS(i,:) = 0.0_wp
            do j = 1,n_atom_gas
               rcm_GAS_uc(i,:) = rcm_GAS_uc(i,:) + r_GAS_uc(i,j,:)*amass(gas_atom(i,j))
               vcm_GAS(i,:) = rcm_GAS_uc(i,:) + GAS_vxyz(i,j,:)*amass(gas_atom(i,j))
               mtmp = mtmp + amass(gas_atom(i,j))
            end do
            rcm_GAS_uc(i,:) = rcm_GAS_uc(i,:)/mtmp
            vcm_GAS(i,:) = vcm_GAS(i,:)/mtmp
         end do
         call msd_samp(n_GAS,rcm_GAS_uc)
         call vcf_samp(n_GAS,vcm_GAS)
         EK = 0.0_wp
         do m = 1, natom
            EK = EK + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         EK_gas = 0.0_wp
         do m = 1, n_GAS
            do i = 1,n_atom_gas
               EK_gas = EK_gas + (0.5_wp/GAS_massi(m,i,1))*dot_product(GAS_vxyz(m,i,:),GAS_vxyz(m,i,:))
            end do
         end do
!print *,'GAS_vxyz 1',GAS_vxyz(1,1,:)
!print *,'GAS_vxyz 1',GAS_vxyz(1,2,:)
!print *,'T_kinetic = ',(EK + EK_gas)*2.0_wp/(K_ev*Ndf)

         write(17,'(8g18.8)') time, (EK + EK_gas), energy + repulenergy*0.5_wp,(EK + EK_gas) + energy + repulenergy*0.5_wp
         !write(18,'(8g18.8)') time, (EK_gas)*2.0_wp/(K_ev*(3*(n_atom_gas*n_GAS) - 3)) ,EK*2.0_wp/(K_ev*(3*natom))
         write(18,'(8g18.8)') time, EK*2.0_wp/(K_ev*(3*natom - 3))

         call flush(17);call flush(18)
         write(*,*) istep,tp(2) - tp(1),tp(3) - tp(2),tp(3) - tp(0)
         if (make_movie) then
            iframe = iframe + 1
            call new_frame_file(imve,'frame',iframe)
            call write_box_frame(imve,n_GAS,print_all = (mod(istep,nprint_all) == 0))
            call new_frame_file(imve,'frame_vel',iframe)
            call write_vel(imve,n_GAS,print_all = (mod(istep,nprint_all) == 0))
         end if
      end if
      end if

!=====================================================end main loop ===
      end do main_loop

      call  msd_print(n_GAS,nsamp,dtreal)
      call  vcf_print(n_GAS,nsamp,dtreal)
      do i = 1,10
         print *,i,sumtp(i)
      end do
      close(17,status='keep')
      close(18,status='keep')
      call write_xyz(24,0)
      close(24,status='KEEP')
CONTAINS

   SUBROUTINE tetra_coords5(bondl,tetra)
      real(wp),intent(in):: bondl
      real(wp),intent(out):: tetra(3,5)
!
! returns a set of coordinates for a tetrahedral molecule
! with the central atom at the origin (atom 1), one atom
! on the z axis below it (atom 2), and atoms 3,4, and 5 on
! a plane above the origin.
!
      real(wp):: zz,sth,cth,rxy,xx,yy
      zz = (1.0_wp/3.0_wp)*bondl
      sth = sqrt(3.0_wp)*0.5_wp
      cth = -0.5_wp
      rxy = (2.0_wp*zz)*sqrt(2.0_wp)
      xx = 0.0_wp
      yy = -rxy
      tetra(1:3,1) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
      tetra(1:3,2) = (/ 0.0_wp, 0.0_wp, -bondl /)
      tetra(1:3,3) = (/ xx, yy, zz /)
      tetra(1:3,4) = (/ (cth*xx - sth*yy), ( sth*xx + cth*yy), zz /)
      tetra(1:3,5) = (/ (cth*xx + sth*yy), (-sth*xx + cth*yy), zz /)
   END SUBROUTINE

END PROGRAM simbox

