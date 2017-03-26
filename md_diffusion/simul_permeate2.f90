
include 'commandline_mod.f90'
include 'precision_mod.f90'
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
include 'msd_mod.f90'
include 'vcf_mod.f90'

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
      integer,parameter:: natom_unit=24
      integer,parameter:: nunitc = 4
      integer,parameter:: number_moves = 251000
      integer,parameter:: n_short = 8
      integer,parameter:: nequib = 1000
      integer,parameter:: nsamp = 25
      integer,parameter:: nprint_all = 2500
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: mN = 14.00_wp
      real(wp),parameter:: mC = 12.00_wp
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,dtreal,qj,sumqj,Temp_K
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp)::CL,time,Ek,rad_Ox(3),rad_N2(3),rad_CO2(3),rr(3),ULJEL,dt_short
      real(wp),allocatable:: massi(:,:),Ox_massi(:,:,:),N2_massi(:,:,:),CO2_massi(:,:,:)
      integer:: j,i,m,ii,ic = 0,natomc,istep,len,status,nat
      integer:: itmp,i_short,narg
      real(wp),allocatable:: coprxyz(:,:)
      real:: timep(0:10),sumtimep(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5
      character(len=132):: infile,carg
      logical:: restart
!
      ntau_msd = 5000
      ntau_vcf = 5000
!
      call cpu_time(timep(0))
      open(unit=24,file='final.pdb')
      open(unit=17,file='energy.out')
      open(unit=1001,file='msd_O2.out')
      open(unit=1002,file='msd_N2.out')
      open(unit=1003,file='msd_CO2.out')
      open(unit=1004,file='vcf_O2.out')
      open(unit=1005,file='vcf_N2.out')
      open(unit=1006,file='vcf_CO2.out')

!
      narg = command_argument_count()
      call get_command_argument (0, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status
         stop
      end if
      write (*,*) 'command name = ', carg(1:len)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) carg(1:len),'  natom  atom_file  T[K]'
         stop
      end if
!
      call get_command_argument (1, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 1
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) natom_max
      print *,'natom = ',natom_max
!
      call get_command_argument (2, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 2
         stop
      end if
      write (*,*) 'arg = ', infile(1:len)
      open(unit=14,file=trim(infile(1:len)))
!
      call get_command_argument (3, carg, len, status)
      if (status /= 0) then
         write (*,*) 'Getting command name failed with status = ', status, ' arg = ', 3
         stop
      end if
      write (*,*) 'arg = ', carg(1:len)
      read(unit=carg,fmt=*) Temp_K
      call get_command_argument (4, infile, len, status)
      if (status /= 0) then
         write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', 4
         restart = .false.
      else
         restart = .true.
         write (*,*) 'restarting from ', infile(1:len)
         open(unit=34,file=trim(infile(1:len)))
      end if

!
!     CL = 0.7168_wp*0.162_wp/0.1552_wp
!     boxl = nunitc*CL
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
!rcutel = boxl2
!rcutel2 = boxl2**2
      rad_Ox = (/sig_O_O2,sig_O_O2,0.0_wp/)*0.5_wp
      rad_N2 = (/sig_N_N2,sig_N_N2,0.0_wp/)*0.5_wp
      rad_CO2 = (/sig_C_CO2,sig_O_CO2,sig_O_CO2/)*0.5_wp
      n_O2 = 40
      n_N2 = 40
      n_CO2 = 40
      natom = natom_max
      T(1) = Temp_K
!
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
      nat = 2
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(Ox_xyz(1:n_O2,1:nat,1:3),N2_xyz(1:n_N2,1:nat,1:3),CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(r_O2_uc(n_O2,nat,3),r_N2_uc(n_n2,nat,3),r_CO2_uc(n_co2,3,3))
      allocate(dr_O2(n_O2,nat,3),dr_N2(n_n2,nat,3),dr_CO2(n_co2,3,3))
      allocate(rcm_O2_uc(n_O2,3),rcm_N2_uc(n_N2,3),rcm_CO2_uc(n_CO2,3))
      allocate(vcm_O2(n_O2,3),vcm_N2(n_N2,3),vcm_CO2(n_CO2,3))
      allocate(Imige_Ox_xyz(1:n_O2,1:3,1:3),Imige_N2_xyz(1:n_N2,1:3,1:3),Imige_CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(fxyz(1:natom_max,3),Ox_fxyz(1:n_O2,1:nat,1:3),N2_fxyz(1:n_N2,1:nat,1:3),CO2_fxyz(1:n_CO2,1:3,1:3))
      allocate(vxyz(1:natom_max,3),Ox_vxyz(1:n_O2,1:nat,1:3),N2_vxyz(1:n_N2,1:nat,1:3),CO2_vxyz(1:n_CO2,1:3,1:3))
      allocate(atom(1:natom_max),Ox_atom(1:n_O2,3),N2_atom(1:n_N2,3),CO2_atom(1:n_CO2,3))
      allocate(charge(natom),Ox_gas_charge(n_O2,3),N2_gas_charge(n_N2,3),CO2_gas_charge(n_CO2,3))
      allocate(proximity(natom_max,4))
      allocate(massi(natom_max,3),Ox_massi(1:n_O2,nat,3),N2_massi(1:n_N2,nat,3),CO2_massi(1:n_CO2,3,3))
      allocate(mol_vec(3,3))
      proximity = 0

      do i = 1,natom
         read (14,'(a6)',advance = 'no') c6
         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
!print *,i,atom(i),'   ',atom_name(atom(i))
         if (atom(i) == iSilicon) massi(i,:) = 1.0_wp/mSi
         if (atom(i) == iOxygen)  massi(i,:) = 1.0_wp/mOx
         if (atom(i) == iOxygenH)  massi(i,:) = 1.0_wp/mOH
      end do

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
! do i = 0,4
! write(888,*) i,qsi(i)
! end do
! do i = 1,natom
! write(888,*) i,atom_name(atom(i)),charge(i),noh(i)
! if (atom(i) ==iSilicon) then
!   sumqj = 0.0_wp
!   do j =1,4
!      ii = proximity(i,j)
!      if (ii == 0) STOP 'error in proximity'
!      qj = charge(ii)
!      if (atom(ii) ==iOxygen) qj = qj*0.5_wp
!      sumqj = sumqj + qj
!   end do
!   write(889,*) i,charge(i),sumqj,charge(i)+sumqj
! end if
! end do
!stop

      call INIT_NLIST(boxl,boxl,boxl,2.61_wp/angstrom)
      call NEW_NLIST

!================== insertion procces==================================
      mol_vec(1,1) =  AOO*0.5_wp
      mol_vec(1,2) = -AOO*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_N2_O2(n_O2,Ox_xyz,1000000,rad_Ox)

      mol_vec(1,1) =  ANN*0.5_wp
      mol_vec(1,2) = -ANN*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_N2_O2(n_N2,N2_xyz,1000000,rad_N2)

      mol_vec(1,1) = 0.0_wp
      mol_vec(1,2) =  ACO2
      mol_vec(1,3) = -ACO2
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_CO2(n_CO2,CO2_xyz,1000000,rad_CO2)
      if (restart) call read_box_frame(34,nmol = 40,print_all = .true.)
! set up the uncorrected coords
      do i = 1, n_O2
         r_O2_uc(i,1,:) = Ox_xyz(i,1,:)
         rr = Ox_xyz(i,2,:) - Ox_xyz(i,1,:)
         call pbc(rr)
         r_O2_uc(i,2,:) = r_O2_uc(i,1,:) + rr
      end do
      do i = 1, n_n2
         r_n2_uc(i,1,:) = N2_xyz(i,1,:)
         rr = N2_xyz(i,2,:) - N2_xyz(i,1,:)
         call pbc(rr)
         r_N2_uc(i,2,:) =  r_n2_uc(i,1,:) + rr
      end do
      do i = 1, n_CO2
         r_CO2_uc(i,1,:) = CO2_xyz(i,1,:)
         rr = CO2_xyz(i,2,:) - CO2_xyz(i,1,:)
         call pbc(rr)
         r_CO2_uc(i,2,:) =  r_CO2_uc(i,1,:) + rr
         !
         rr = CO2_xyz(i,3,:) - CO2_xyz(i,1,:)
         call pbc(rr)
         r_CO2_uc(i,3,:) = r_CO2_uc(i,1,:) + rr
      end do


      if (make_movie) then ! write out first frame of movie
         iframe = 0
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_N2,print_all = .true.)
      end if
!========================================================================
!

      do i = 1,n_O2
        Ox_gas_charge(i,1) = q_O_O2
        Ox_gas_charge(i,2) = q_O_O2
        Ox_gas_charge(i,3) = -2.0_wp*q_O_O2
        Ox_atom(i,1) = iO_O2
        Ox_atom(i,2) = iO_O2
        Ox_massi(i,:,:) = 1.0_wp/mOx
      end do
      do i = 1,n_N2
        N2_gas_charge(i,1) = q_N_N2
        N2_gas_charge(i,2) = q_N_N2
        N2_gas_charge(i,3) = -2.0_wp*q_N_N2
        N2_atom(i,1) = iN_N2
        N2_atom(i,2) = iN_N2
        N2_massi(i,:,:) = 1.0_wp/mN
      end do
      do i = 1,n_CO2
        CO2_gas_charge(i,1) = q_C_CO2
        CO2_gas_charge(i,2) = q_O_CO2
        CO2_gas_charge(i,3) = q_O_CO2
        CO2_atom(i,1) = iC_CO2
        CO2_atom(i,2) = iO_CO2
        CO2_atom(i,3) = iO_CO2
        CO2_massi(i,1,:) = 1.0_wp/mC
        CO2_massi(i,2:3,:) = 1.0_wp/mOx
      end do

      call LJ_INIT

!     call comvel(natom,n_co2,n_o2,n_n2,t(1),vxyz,ox_vxyz,n2_vxyz,co2_vxyz)
!
!     vxyz = vxyz*sqrt(massi)
!     ox_vxyz = ox_vxyz*sqrt(ox_massi)
!     n2_vxyz = n2_vxyz*sqrt(n2_massi)
!     co2_vxyz = co2_vxyz*sqrt(co2_massi)
      vxyz = 0.0_wp
      Ox_vxyz = 0.0_wp
      N2_vxyz = 0.0_wp
      CO2_vxyz = 0.0_wp

      call set_bond_angle_lists()
      call msd_init(40)
      call vcf_init(40)

      time = 0.0_wp
      dtreal = 2.0_wp
      dt = dtreal*dt_fs
      dt_short = dt/n_short

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      energy = 0.0_wp
!
!      call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!      call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!      call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! Electrostatic
!     call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)
!     call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_n2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)
!     call FORCE_ENERGY_EL(1,natom,1,n_co2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)

!======================================================== main loop ===
      main_loop:do istep = 1,number_moves

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi
      Ox_vxyz = Ox_vxyz + 0.5_wp*dt*Ox_fxyz*Ox_massi
      N2_vxyz = N2_vxyz + 0.5_wp*dt*N2_fxyz*N2_massi
      CO2_vxyz = CO2_vxyz + 0.5_wp*dt*CO2_fxyz*CO2_massi

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      call cpu_time(timep(1))
      call force_keating
      call repul_energy2(1,natom)
      call O2_stretching
      call N2_stretching
      call CO2_stretching
      call CO2_angle_bending

!=============================== short loop ===
      short_loop: do i_short = 1, n_short
         coprxyz = rxyz
         call vverlet_a(dt_short,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)

         call gas_vverlet_a(dt_short,Ox_massi,Ox_xyz,Ox_vxyz,Ox_fxyz,dr_O2)
         call pbc(Ox_xyz)
         r_O2_uc = r_O2_uc + dr_O2

         call gas_vverlet_a(dt_short,N2_massi,N2_xyz,N2_vxyz,N2_fxyz,dr_N2)
         call pbc(N2_xyz)
         r_N2_uc = r_N2_uc + dr_N2

         call gas_vverlet_a(dt_short,CO2_massi,CO2_xyz,CO2_vxyz,CO2_fxyz,dr_CO2)
         call pbc(CO2_xyz)
         r_CO2_uc = r_CO2_uc + dr_CO2

         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3))) then
               !! print *,'NEW LIST'
               call NEW_NLIST
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         Ox_fxyz = 0.0_wp
         N2_fxyz = 0.0_wp
         CO2_fxyz = 0.0_wp
         energy = 0.0_wp
!call cpu_time(timep(4))
         call force_keating
!call cpu_time(timep(5)); sumtimep(4) = sumtimep(4) + timep(5)-timep(4)
         call repul_energy2(1,natom)
!call cpu_time(timep(6)); sumtimep(5) = sumtimep(5) + timep(6)-timep(5)
         call O2_stretching
         call N2_stretching
         call CO2_stretching
!call cpu_time(timep(7)); sumtimep(6) = sumtimep(6) + timep(7)-timep(6)
         call CO2_angle_bending
!call cpu_time(timep(8)); sumtimep(7) = sumtimep(7) + timep(8)-timep(7)

         call vverlet_b(dt_short,massi,vxyz,fxyz)
         call gas_vverlet_b(dt_short,Ox_massi,Ox_vxyz,Ox_fxyz)
         call gas_vverlet_b(dt_short,N2_massi,N2_vxyz,N2_fxyz)
         call gas_vverlet_b(dt_short,CO2_massi,CO2_vxyz,CO2_fxyz)
         !coprxyz = rxyz
!call cpu_time(timep(9)); sumtimep(8) = sumtimep(8) + timep(9)-timep(8)
!============================end short loop ===
      end do short_loop
      call cpu_time(timep(2)); sumtimep(1) = sumtimep(1) + timep(2) - timep(1)

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
!
!     call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!     call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!     call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL)
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! Electrostatic
!     call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)
!     call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_n2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)
!     call FORCE_ENERGY_EL(1,natom,1,n_co2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)

      call cpu_time(timep(3)); sumtimep(2) = sumtimep(2) + timep(3) - timep(2)

      vxyz = vxyz + 0.5_wp*dt*fxyz*massi
      Ox_vxyz = Ox_vxyz + 0.5_wp*dt*Ox_fxyz*Ox_massi
      N2_vxyz = N2_vxyz + 0.5_wp*dt*N2_fxyz*N2_massi
      CO2_vxyz = CO2_vxyz + 0.5_wp*dt*CO2_fxyz*CO2_massi

      time = time + dtreal
      if (istep <= nequib) then
         !  print *,sum(fxyz(1:natom,:),dim=1)
         Ek = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         do m = 1, n_O2
            do i = 1,2
               Ek = Ek + (0.5_wp/Ox_massi(m,i,1))*dot_product(Ox_vxyz(m,i,:),Ox_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_N2
            do i = 1,2
               Ek = Ek + (0.5_wp/N2_massi(m,i,1))*dot_product(N2_vxyz(m,i,:),N2_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_CO2
            do i = 1,3
               Ek = Ek + (0.5_wp/CO2_massi(m,i,1))*dot_product(CO2_vxyz(m,i,:),CO2_vxyz(m,i,:))
            end do
         end do

         T(2) = Ek*2.0_wp/( 3.0_wp*K_ev*(natom + 280) )
         if (mod(istep,nsamp) == 0) then
            write(*,*) istep,timep(2) - timep(1),timep(3) - timep(2),timep(3) - timep(0)
            write(17,'(4g18.8)') time, EK, energy + repulenergy*0.5_wp, Ek + energy + repulenergy*0.5_wp
            call flush(17)
         end if
         vxyz = vxyz*sqrt(T(1)/T(2))
         Ox_vxyz = Ox_vxyz*sqrt(T(1)/T(2))
         N2_vxyz = N2_vxyz*sqrt(T(1)/T(2))
         CO2_vxyz = CO2_vxyz*sqrt(T(1)/T(2))
      end if

      if ((mod(istep,nsamp) == 0) .AND. (istep >= nequib)) then
         forall(i = 1:n_O2)
            rcm_O2_uc(i,:) = ( r_O2_uc(i,1,:) + r_O2_uc(i,2,:) )*0.5_wp
            rcm_N2_uc(i,:) = ( r_N2_uc(i,1,:) + r_N2_uc(i,2,:) )*0.5_wp
            rcm_CO2_uc(i,:) = ( r_CO2_uc(i,1,:)*mC + r_CO2_uc(i,2,:)*mOx + r_CO2_uc(i,3,:)*mOx )/(mOx + mOx + mC)
            vcm_O2(i,:) = ( Ox_vxyz(i,1,:) + Ox_vxyz(i,2,:) )*0.5_wp
            vcm_N2(i,:) = ( N2_vxyz(i,1,:) + N2_vxyz(i,2,:) )*0.5_wp
            vcm_CO2(i,:) = ( CO2_vxyz(i,1,:)*mC + CO2_vxyz(i,2,:)*mOx + CO2_vxyz(i,3,:)*mOx )/(mOx + mOx + mC)
         end forall
         call msd_samp(n_O2,rcm_O2_uc,rcm_N2_uc,rcm_CO2_uc)
         call vcf_samp(n_O2,vcm_O2,vcm_N2,vcm_CO2)
         Ek = 0.0_wp
         do m = 1, natom
            Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
         end do
         do m = 1, n_O2
            do i = 1,2
               Ek = Ek + (0.5_wp/Ox_massi(m,i,1))*dot_product(Ox_vxyz(m,i,:),Ox_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_N2
            do i = 1,2
               Ek = Ek + (0.5_wp/N2_massi(m,i,1))*dot_product(N2_vxyz(m,i,:),N2_vxyz(m,i,:))
            end do
         end do
         do m = 1, n_CO2
            do i = 1,3
               Ek = Ek + (0.5_wp/CO2_massi(m,i,1))*dot_product(CO2_vxyz(m,i,:),CO2_vxyz(m,i,:))
            end do
         end do

         write (17,'(4g18.8)') time, EK, energy + repulenergy*0.5_wp, Ek + energy + repulenergy*0.5_wp
         call flush(17)
         write(*,*) istep,timep(2) - timep(1),timep(3) - timep(2),timep(3) - timep(0)
         if (make_movie) then
            iframe = iframe + 1
            call new_frame_file(imve,'frame',iframe)
            call write_box_frame(imve,n_n2,print_all = (mod(istep,nprint_all) == 0))
         end if
      end if

!=====================================================end main loop ===
      end do main_loop
      call  msd_print(40,nsamp,dtreal)
      call  vcf_print(40,nsamp,dtreal)
      do i = 1,10
         print *,i,sumtimep(i)
      end do
      close(17,status='keep')
      close(1001,status='keep')
      close(1002,status='keep')
      close(1003,status='keep')
      close(1004,status='keep')
      close(1005,status='keep')
      close(1006,status='keep')
      call write_xyz(24,0)
      close(24,status='KEEP')

END PROGRAM simbox

