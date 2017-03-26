
INCLUDE 'precision_mod.f90'
INCLUDE 'rand_mod.f90'
INCLUDE 'M_intcat_jlg.f90'
INCLUDE 'misc_legacy_mod.f90'
INCLUDE 'math_const_mod.f90'
INCLUDE 'phys_const_mod.f90'
INCLUDE 'command_line_mod.f90'
INCLUDE 'command_arg_mod.f90'
INCLUDE 'files_mod.f90'
INCLUDE 'sort_mod.f90'
INCLUDE 'global_vars_mod_KMC.f90'
INCLUDE 'grid_data_mod.f90'
!INCLUDE 'netcdf_array_mod.f90'
INCLUDE 'list_mod.f90'
INCLUDE 'listset_mod.f90'
INCLUDE 'sets.f90'
INCLUDE 'HKNonLattice2.f90'
INCLUDE 'hk_lattice_mod.f90'
!INCLUDE 'hk_lattice_nonlattice_mod.f90'
INCLUDE 'hk_lattice_nonlattice2_mod.f90'
INCLUDE 'unionfind_mod.f90'
INCLUDE 'hk_lattice_uf_mod.f90'
INCLUDE 'percolation_lattice_mod.f90'
INCLUDE 'Keating_parameters_mod_vonAlfthan.f90'
!INCLUDE 'Dreiding_parameters_mod.f90'
!INCLUDE 'charmm_mod.f90'
INCLUDE 'seaton_mod.f90'
INCLUDE 'atom_types_mod.f90'
!INCLUDE 'atom_types_mod_TEOS.f90'
INCLUDE 'constants_mod.f90'
!INCLUDE 'constants_mod_TEOS.f90'
INCLUDE 'coordinates_mod.f90'
!INCLUDE 'coordinates_mod_2DP.f90'
INCLUDE 'connectivity_mod.f90'
INCLUDE 'check_structure_mod_Silica.f90'
!INCLUDE 'check_structure_mod_TEOS.f90'
INCLUDE 'teos_atomlst_mod.f90'
INCLUDE 'nlist_mod3.f90'
INCLUDE 'verlet_list_nl.f90'
INCLUDE 'find_atoms_mod.f90'
INCLUDE 'charges2_mod.f90'
!INCLUDE 'charges2_mod_TEOS.f90'
INCLUDE 'matvec3d_mod.f90'
INCLUDE 'quat2mat_mod.f90'
INCLUDE 'rotate_axis_mod.f90'
INCLUDE 'tetra_coords_mod.f90'
INCLUDE 'probe_mol_mod.f90'
INCLUDE 'probe_mol_init_mod.f90'
!INCLUDE 'probe_mol_init_mod_TEOS.f90'
INCLUDE 'bond_angle_types_mod.f90'
!INCLUDE 'bond_angle_types_mod_TEOS.f90'
INCLUDE 'bond_angle_list_mod.f90'
!INCLUDE 'box_frames_mod.f90'
INCLUDE 'frames_mod.f90'
!INCLUDE 'comvel_mod.f90'
!INCLUDE 'vverlet_mod.f90'
!INCLUDE 'glassbond_mod.f90'
INCLUDE 'repul_force_NL_vonAlfthan.f90'
INCLUDE 'ran_point_sphere_mod.f90'
INCLUDE 'insert_mod_2.f90'
INCLUDE 'lj_el_mod_O2_N2_CO2.f90'
INCLUDE 'force_en_el_mod.f90'
INCLUDE 'bonded_force_gas_mod_O2_N2_CO2.f90'
!INCLUDE 'tcf_mod.f90'
!INCLUDE 'msd_mod3.f90'
!INCLUDE 'vcf_mod3.f90'
INCLUDE 'rdf_mod.f90'
INCLUDE 'utility_mod.f90'
INCLUDE 'energy_keating_mod.f90'
INCLUDE 'force_keating_mod.f90'
INCLUDE 'statistics_mod.f90'
INCLUDE 'timer_mod.f90'
INCLUDE 'force_energy_minimize_mod.f90'
INCLUDE 'Steepest_Descent_mod.f90'
INCLUDE 'relax_mc_mod.f90'
INCLUDE 'voidage5_mod.f90'
INCLUDE 'Henrys_law_calc_mod.f90'
INCLUDE 'density_profile_mod.f90'
INCLUDE 'deposition5_mod.f90'
INCLUDE 'relax_VL_mod.f90'
INCLUDE 'select_atom_mod.f90'
INCLUDE 'attach_tetrahedra_mod.f90'
!INCLUDE 'attach_molecule_mod.f90'
!INCLUDE 'attach_molecule_mod_TEOS.f90'
INCLUDE 'water_attack_VL_mod.f90'
!INCLUDE 'water_attack_VL_mod_H.f90'
INCLUDE 'condensation_mod.f90'
!INCLUDE 'condensation_mod_H.f90'
INCLUDE 'surface_atoms_mod.f90'
INCLUDE 'switch_mod.f90'
INCLUDE 'relax_bonds_VL_NL_mod.f90'
INCLUDE 'list_atoms_in_cells_mod.f90'
!INCLUDE 'list_atoms_in_cells_mod_TEOS.f90'


PROGRAM porous_silica
      USE precision_mod
      USE command_line_mod
      USE command_arg_mod
      USE phys_const_mod
      USE global_vars_mod
      USE constants_mod
      USE atom_types_mod
      USE list_mod
      USE connectivity_mod
      USE check_structure_mod
      USE coordinates_mod
      USE rand_mod
      USE frames_mod
      USE Relax_mod
      USE nlist_mod
      USE verlet_list_mod
      USE HKNonLattice_mod
      USE density_profile_mod
      USE deposition_mod
      USE rdf_mod
      USE Keating_parameters_mod, only: ASiO
      USE water_attack_mod
      USE files_mod
      USE Relax_mod
      USE select_atom_mod
      USE attach_molecule_mod
      USE probe_mol_mod
      USE probe_mol_init_mod
      USE voidage_mod
      USE surface_atoms_mod
      USE energy_keating_mod
      USE repul_energy_mod
      USE seaton_mod
      USE charges_mod
      USE lj_el_mod
      USE Henrys_law_calc_mod
      USE condensation_mod
      USE relax_bonds_mod
      USE find_atoms_mod
      USE list_atoms_in_cells_mod
      USE ran_point_sphere
      USE quat2mat_mod
      USE M_intcat_jlg
      USE statistics_mod
      implicit none
      logical:: make_movie
      real(wp):: top_atom,void_crit,rverlet,dlayer,Hmax
      real(wp):: dLflexible, dLfixed, dlkeep,quat(4),aa(3,3),ulj,uel,r1,rc,duel
      real(wp):: Rbox_min(3),Rbox_max(3),del_nlist, fover,rl(3),ru(3),energy
      real(wp):: dens_gcm3,volbox,rcutH,dl,Uljel_1,Khenry,facc,qj,sumqj,rv(3)
      integer:: i,itmp,j,ipconfig,ntotal,imve,narg,nc,ii,iz,k,it
      real:: timep(0:10)
      character:: date*8,ctime*10
      character(len=32):: ctmp,c6*6,c5*5,c32
      character(len=132):: coordfile,restart_coordfile,carg,ftyp*3,ctrlfile
      integer:: io_config, io_config_final, io_topol, io_ctrl, io_restart_xyz
      logical:: success,ok,booltmp,restart
      real(wp),allocatable:: massi(:,:),mass(:,:)
      integer,allocatable:: listOH(:),listALL(:),listSph(:)
      integer:: nSph,nlst,nfixed,noxdel,n_oh,n_si,n_ob
      type(probe_mol):: pr,pr1
      type(statistic):: s_lj,s_el,s_tot
integer,parameter:: nrxn = 4
integer,parameter:: nbinmax = 10000
real(wp):: rate(nrxn,nbinmax),rate_sum,psum
integer:: irxn(nrxn,nbinmax,0:3)
type(mylist):: alst(3),oh_lst,ohb_lst,lst0,lstA
integer:: ia_OH
integer:: ia_OSiOH3
integer:: ia_OSiO
integer:: ib,ir,iOH,stat
real(wp):: time,delta_time,xi1
real(wp),parameter:: rohbl = 4.0_wp/Angstrom
real(wp),parameter:: volohbl = (4.0_wp/3.0_wp)*pi*(rohbl**3)
integer:: n_neigh,neighbors(2000)
real(wp):: P_SiOH4,P_H2O,exp_Ea_kt,T_Kelvin,zu,zl,vco2,vo2,vn2
logical:: hydrothermal = .false.
!
      print *,'Kinetic Monte Carlo'
      print *,'Number of dimensions = ',NDIM
      print *,'Number of periodic dimensions = ',NPER
      if (NDIM < NPER) then
         print *,'The number of dimensions must be >= the number of periodic dimensions'
         STOP 'Change NDIM &/ NPER in coordinates_mod and recompile'
      end if
      call cpu_time(timep(0))
!-----------------------------------------------------------
! Read in simulation data
!     kmc_step = number of kmc steps so far
!     ntotal = total number of kmc steps
!     irand_seed = random number seed
!-----------------------------------------------------------
      narg = command_argument_count()
      call next_command_argument(carg, "command",stat)
      if (narg < 6) then
         write(*,*)'usage :'
         write(*,*) trim(carg),'control_file  atom_file  kmc_step  ntotal  T[K]  P_H2O/P_Si(OH)4  random_seed'
         stop
      end if
      call next_command_argument(ctrlfile, "control file", stat)
      call next_command_argument(coordfile, "input file",stat)
      nc = len_trim(coordfile)
      ftyp = coordfile(nc - 2:nc)
      call get_free_file_unit(io_config)
      open(unit=io_config,file=trim(coordfile))

      call next_command_argument(kmc_step, "KMC step", stat)
      call next_command_argument(ntotal, "ntotal", stat)
      call next_command_argument(T_kelvin, "Temperature [K]", stat)
      call next_command_argument(fwater, "fwater", stat)
      call next_command_argument(irand_seed, "irand_seed", stat)

      call get_free_file_unit(io_ctrl)
      open(unit=io_ctrl,file=trim(ctrlfile))
!      read(io_ctrl,*) c32, kmc_step            ; write(*,*) trim(c32), kmc_step 
!      read(io_ctrl,*) c32, ntotal              ; write(*,*) trim(c32), ntotal
!      read(io_ctrl,*) c32, T_kelvin            ; write(*,*) trim(c32), T_kelvin
!      read(io_ctrl,*) c32, fwater              ; write(*,*) trim(c32), fwater
!      read(io_ctrl,*) c32, irand_seed          ; write(*,*) trim(c32), irand_seed
      read(io_ctrl,*) c32,void_crit            ; write(*,*) trim(c32), void_crit
      read(io_ctrl,*) c32,P_SiOH4              ; write(*,*) trim(c32), P_SiOH4
      read(io_ctrl,*) c32, del_nlist           ; write(*,*) trim(c32), del_nlist
      read(io_ctrl,*) c32, rverlet             ; write(*,*) trim(c32), rverlet
      read(io_ctrl,*) c32, NRELAX              ; write(*,*) trim(c32), NRELAX
      read(io_ctrl,*) c32, Hmax                ; write(*,*) trim(c32), Hmax
      read(io_ctrl,*) c32, e_activ             ; write(*,*) trim(c32), e_activ
      read(io_ctrl,*) c32, del_rxn             ; write(*,*) trim(c32), del_rxn
      read(io_ctrl,*) c32, ipconfig            ; write(*,*) trim(c32), ipconfig
      read(io_ctrl,*) c32, make_movie          ; write(*,*) trim(c32), make_movie
      read(io_ctrl,*) c32, Hydrothermal        ; write(*,*) trim(c32), Hydrothermal
!      read(io_ctrl,*) c32, fover               ; write(*,*) trim(c32), fover
      del_nlist = del_nlist/Angstrom
      rverlet = rverlet/Angstrom
      Hmax = Hmax/angstrom
      close(io_ctrl)
!
!      del_nlist = 0.5
!      P_SiOH4 = 1333.22375
!      void_crit = 0.035
!      ipconfig = 50
!      Hmax = 150.0_wp/Angstrom
!      NRELAX = 50
!      rverlet = 0.5
!      e_activ = 0.75
!      del_rxn = 0.14
!      Hydrothermal = .false.
!
      kBoltzT = K_ev*T_Kelvin
      exp_Ea_kt = exp(-del_rxn/kBoltzT)
      P_SiOH4 = P_SiOH4/upressure
      P_H2O = fwater*P_SiOH4
      if(hydrothermal) then
         P_H2O = P_SiOH4
         P_SiOH4 = 0.0_wp
      end if
!
      print '(a,g18.9,a)','ulength   = ',ulength,' m'
      print '(a,g18.9,a)','uenergy   = ',uenergy,' J'
      print '(a,g20.10,a)','upressure = ',upressure,' Pa'
      print '(a,g18.9,a,g18.9,a)','P_SiOH4 = ',P_SiOH4,' = ',P_SiOH4*upressure,' Pa'
      print '(a,g18.9,a)',        '        = ',P_SiOH4*upressure/torr2Pa,' torr'
      print '(a,g18.9,a)','kT   = ',kBoltzT,' eV'
      print *,'T = ',T_Kelvin,' K'
      print '(a,g18.9)','exp(-EA/kT) = ',exp_Ea_kt
!
! Read in the initial molecular system file header
      select case(ftyp)
      case('xyz')
         call read_xyz_header(io_config,natom,boxl)
         io_topol = 15
         call get_free_file_unit(io_topol)
         open(unit=io_topol,file=trim(coordfile(1:nc - 3))//'con')
      case('pdb')
         call read_pdb_header(io_config,natom,boxl)
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
!
      natom_max = natom
      !natom_max = natom + 7*ntotal  ! 3d
      call SET_BOX_SIZE(boxl)
      volbox = boxl(1)*boxl(2)*boxl(3)*(Angstrom**3) ! Box volume in A^3
!
! Allocate storage for solid
      allocate(rxyz(natom_max,3))
      allocate(fxyz(natom_max,3))
      allocate(atom(natom_max))
      atom = 0
!     allocate(atom(natom_max),latom(natom_max))
!     atom2 = 0; latom = .false.
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
      proximity = 0

!     allocate(coprxyz(natom_max,3))
!     allocate(vxyz(1:natom_max,3))
!
      allocate(mass(natom_max,3),massi(natom_max,3))

! read in coordinates and atom types
      select case(ftyp)
      case('xyz')
         call read_xyz(io_config,natom,atom,rxyz)
         close(io_config)
      case('pdb')
         call read_pdb(io_config,natom,atom,rxyz)
      end select

! read in connectivity
      select case(ftyp)
      case('xyz')
         call read_conn(io_topol,proximity)
         close(io_topol)
      case('pdb')
         call read_conect_pdb(io_config,natom,proximity)
         close(io_config)
      end select

! assign mases and charges
      do i = 1,natom
         mass(i,:) = amass(atom(i))
      end do
      massi(1:natom,:) = 1.0_wp/mass(1:natom,:)
      call assign_charge(1,natom)

      write(*,'(a,f12.8)') 'sum charge = ',sum(charge(1:natom))
      write(*,'(a,f12.8)') 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)

      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'N_Si = ',n_si/volbox,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob/volbox,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh/volbox,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*(1e24_wp/(volbox*NAvo))
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'

      call LJ_init
      call init_probe_mols()

!   allocate(listOH(count(atom == iOxygenH)))
!   n_OH = 0
!   do i = 1,natom
!      if (atom(i) == iOxygenH) then
!         n_OH = n_OH + 1
!         listOH(n_OH) = i
!      end if
!   end do
!   print *,'n_OH = ',n_OH

! check the structure
      call check_proximity(1,natom,ok,i)
      if (.NOT.ok) then
         write(*,*) 'proximity array is NOT consistent'
         write(*,*) 'check_proximity ',ok,i
         call write_atom_info(6,i)
         stop
      end if
      call check_proximity2(1,natom,atom,proximity,ok,i)
      if (.NOT.ok) then
        write(*,*) 'proximity array is NOT consistent'
        write(*,*) 'check_proximity2 ',ok,i
        call write_atom_info(6,i)
        stop
      end if

! Initially the box should be centered at the origin.
      print '(a,3f12.6)','boxl old = ',boxl(3)

      ! translate the box by -L_z/2
      rxyz(1:natom,3) = rxyz(1:natom,3) - boxl2(3)

      ! and set the new Lz to Max(Hmax,2*Lz)
      boxl(3) = max(2.0_wp*boxl(3),Hmax)
      call SET_BOX_SIZE(boxl)
      print '(a,3f12.6)','    boxl = ',boxl(3)
      print '(a,3f12.6)','   boxl2 = ',boxl2(3)


!     call new_frame_file(imve,'frame',10000,'pdb')
!     call write_frame(imve,1,natom)
      call new_frame_file(imve,'frame',10000,'xyz')
      call write_xyz_file(imve,1,natom,atom,rxyz,boxl,box_corners=.true.,comment=8//' Box corners added')
      call new_frame_file(imve,'frame',10000,'con')
      call write_conn(imve,1,natom,proximity,'connectivity at step = '//j)


      dLflexible = 5.0_wp/Angstrom
      dLfixed = 5.0_wp/Angstrom
      dlkeep = dLflexible + dLfixed

      Rbox_min = -boxl2
      Rbox_min(3)=Rbox_min(3)-5.0_wp/Angstrom
      !Rbox_min = (/ -boxl2(1), -boxl2(2), -dLflexible /)
      Rbox_max = (/  boxl2(1),  boxl2(2), Hmax /)

      !call INIT_NLIST(Rbox_min,Rbox_max,del_nlist,natom_max)
      call INIT_NLIST(Rbox_min,Rbox_max,del_nlist,natom_max)
      if (del_nlist < maxval(sigLJ)) STOP 'Neighbor list is too small'
      call Init_Verlet_List(rv0=rverlet,rcut0=2.6_wp/Angstrom)
      call NEW_NLIST(1,natom)
      call new_vlist

boxl=1000.0_wp
boxl2 = boxl/2.0_wp
boxli = 1.0_wp/boxl

rv = (/ 0.0_wp, 0.0_wp, -2.5_wp /)/Angstrom
do iz = 1,10
rv(3) = rv(3) + 5.0_wp/Angstrom
call reset_stat(s_lj)
call reset_stat(s_el)
call reset_stat(s_tot)
do it = 1,200
pr = probe_sioh4
pr1 = pr
!call random_rotation(aa)
call ran4sph(quat)
call quat2mat(quat,aa)
pr1%r = matmul(aa,pr%r)
pr1%r(1,:) = pr1%r(1,:) + rv(1)
pr1%r(2,:) = pr1%r(2,:) + rv(2)
pr1%r(3,:) = pr1%r(3,:) + rv(3)
do i = 1,pr1%n
   call pbc(pr1%r(1:3,i))
end do
call ENERGY_LJ_EL2(pr1,1,natom,Ulj,Uel)
call update_stat(s_lj,Ulj)
call update_stat(s_el,Uel)
call update_stat(s_tot,Uel+Ulj)

write(88,'(i8,9(1x,f0.8))'),it,rv(3)*Angstrom,ulj,uel,ulj+uel,s_lj%mean,s_el%mean,s_tot%mean
end do
write(88,'(/)')

end do
call new_list(lstA,natom)
call new_list(lst0,(5*natom_max/ncelz))
print *,'ncelz = ',ncelz
do iz = ncelz,1,-1
   call Z_NEIGH_CELL(iz,n_neigh,neighbors)
   call list_atoms_in_cells(n_neigh,neighbors,lst0)
   call add_to_list(lstA,lst0)

   sumqj = 0.0_wp
   do j = 1,lstA%n
      qj = charge(lstA%ind(j))
      sumqj = sumqj + qj
   end do
   call ENERGY_LJEL_LIST_NEUTRAL(pr1,lstA%n,lstA%ind,atom,charge,rxyz,Ulj,Uel)

   rc = -boxl2(3) - (iz-1)*del_nl_z
   duel = 0.0_wp
   do k = 1,pr%n
      r1 = pr%r(3,k) - rc
      duel = duel + pr%q(k)*sumqj/r1
   end do
if (lstA%n==0) cycle
   print '(3i9,6(1x,g16.8))',iz,neighbors(1),lstA%n,sumqj/qi(iOxygen),ulj,uel,duel,uel-duel
end do

stop

      call Init_HKNonLattice(natom_max)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
      write(*,*) 'n_cluster = ',n_cluster
!     call hk_initialize((nb(1)*nb(2)*nb(3))/2)
!
! output files
!      OPEN(unit= 15,file= 'psil_config.pdb')
!      OPEN(unit= 20,file= 'psil_topology.out')
!      OPEN(unit= 21,file= 'psil_rand.out')
!      OPEN(unit= 22,file= 'psil_dens.out')
!      OPEN(unit= 23,file= 'psil_void_SiO4.out')
!      OPEN(unit= 24,file= 'psil_topatom.out')
!      OPEN(unit= 25,file= 'psil_time.out')
!      OPEN(unit= 29,file= 'psil_kmc_history.out')
!      OPEN(unit= 77,file= 'psil_SiO4_attempt.out')
!      OPEN(unit= 78,file= 'psil_SiO4_add.out')
!      OPEN(unit= 79,file= 'psil_H2O_add.out')
!      OPEN(unit= 80,file= 'psil_H2O_remove.out')
!      OPEN(unit= 81,file= 'psil_bond_switch.out')
!
      call cpu_time(timep(0))
!
!----------------------------------------------------------------------
!      time = 0.0_wp
!      main_loop: do while (kmc_step < ntotal)
!      call cpu_time(timep(1))
!      kmc_step = kmc_step + 1
!!
!      call assign_charge(1,natom)
!
!      dlayer = (5.0_wp/angstrom)
!      volbin = boxl(1)*boxl(2)*dlayer
!!      call henry_profile_z(probe_SiO4,dlayer,1)
!!      call henry_profile_z(probe_tip3p,dlayer,2)
!      call henry_profile(probe_SiO4,dlayer,1)
!      call henry_profile(probe_tip3p,dlayer,2)
!      call den_profile_rxn(1,natom,dlayer)
!!      call prob_dist(henryk(1:nbin,1),inthenryk(1:nbin,1),void_crit,itop_nonp_SiO4)
!!      call prob_dist(henryk(1:nbin,2),inthenryk(1:nbin,2),void_crit,itop_nonp_H2O)
!!do i = 1,nbin
!!   write(99,'(i5,2(1x,f0.6))')i,henryk(i,1),henryk(i,2)
!!end do
!!write(99,*)
!!
!!                         k_1
!! (1)  |-OH  +  Si(OH)4   --->   |-Si(OH)3  +  H2O
!!
!!                         k_1R
!! (2)  |-Si(OH)3  +  H2O  --->   |-OH  +  Si(OH)4
!!
!!                         k_2
!! (3)     |-OH  +  HO-|   --->   |-O-|  +  H2O
!!
!!                         k_2R
!! (4)     |-O-|  +  H2O   --->   |-OH  +  HO-|
!!
!      rate_sum = 0.0_wp
!      rate(1:nrxn,1:nbin) = 0.0_wp
!      irxn(1:nrxn,1:nbin,:) = 0
!      do ib = 1,nbin
!         call Z_NEIGH_CELL(ib,n_neigh,neighbors)
!         call list_atoms_in_cells(n_neigh,neighbors,alst)
!         call rand_from_list(alst(L_OH),ia_OH)
!         call rand_from_list(alst(L_OSiOH3),ia_OSiOH3)
!         call rand_from_list(alst(L_OSiO),ia_OSiO)
!         if(ia_OSiO /= 0)then
!            rate(4,ib) = exp_Ea_kt*henryk(ib,2)*P_H2O*ndenrxn(ib,iOSiO)
!            if(ib < itop_nonp_H2O) rate(4,ib)=0.0_wp
!            rate_sum = rate_sum + rate(4,ib)
!            irxn(4,ib,0) = 1
!            irxn(4,ib,1) = ia_OSiO
!         end if
!         if(ia_OSiOH3 /= 0)then
!            rate(2,ib) = exp_Ea_kt*henryk(ib,2)*P_H2O*ndenrxn(ib,iOSiOH3)
!            if(ib < itop_nonp_H2O) rate(2,ib) = 0.0_wp
!            rate_sum = rate_sum + rate(2,ib)
!            irxn(2,ib,0) = 1
!            irxn(2,ib,1) = ia_OSiOH3
!         end if
!         if(ia_OH /= 0)then
!            rate(1,ib) = henryk(ib,1)*ndenrxn(ib,iOxygenH)*P_SiOH4
!            if(ib < itop_nonp_SiO4) rate(1,ib) = 0.0_wp
!            rate_sum = rate_sum + rate(1,ib)
!            irxn(1,ib,0) = 1
!            irxn(1,ib,1) = ia_OH
!            call find_atoms(ia_OH, iOxygenH, rohbl, oh_lst%n, oh_lst%ind)
!            ohb_lst%n = 0
!            ohb_lst%ind = 0
!            do i = 1,oh_lst%n
!               iOH = oh_lst%ind(i)
!               if(OH_groups_OK(ia_OH,iOH)) call add_to_list(ohb_lst,iOH)
!            end do
!            call rand_from_list(ohb_lst,iOH)
!            if(iOH /= 0)then
!               rate(3,ib) = kBoltzT*(ndenrxn(ib,iOxygenH)/volbin)*ohb_lst%n*(volbin/volohbl)
!               if(ib <  itop_nonp_H2O) rate(3,ib) = 0.0_wp
!               rate_sum = rate_sum + rate(3,ib)
!               irxn(3,ib,0) = 2
!               irxn(3,ib,1) = ia_OH
!               irxn(3,ib,2) = iOH
!            end if
!         end if
!      end do
!
!      xi1 = rand()
!      psum = 0.0_wp
!      bin_loop: do ib = 1,nbin
!      do ir = 1,nrxn
!         psum = psum + rate(ir,ib)/rate_sum
!         if (psum >= xi1) then
!            exit bin_loop
!         end if
!      end do
!      end do bin_loop
!      write(*,'(g18.8,i4,2i9)') rate_sum,ir,ib,nbin
!      write(29,'(g18.8,i4,2i9)') rate_sum,ir,ib,nbin
!!
!!     carry out reaction #ir in bin #ib
!print *,'ir = ',ir
!
!do
!   ia_OH = int(rand()*natom) + 1
!   if (atom(ia_OH) == iOxygenH) exit
!end do
!call find_atoms(ia_OH, iOxygenH, rohbl, oh_lst%n, oh_lst%ind)
!call rand_from_list(ohb_lst,iOH)
!if (ir > 4)then
!if(rand()<0.9)then
!    ir = 1
!else
!    ir = 3
!  if(iOH==0)ir=1
!end if
!end if
!print *,'ir = ',ir
!
!      select case(ir)
!      case(1) ! reaction 1
!!         ia_OH = irxn(1,ib,1)
!         ndel = 0
!         idel(1) = 0
!         !call attach_tetrahedra(ia_OH,aSiOSi,aSiO,overlap_check=.false.,success=success)
!         call attach_mol( ia,ndel,idel,bondl_SiO,aSiOSi,4,transpose(probe_SiO3%r), &
!                             probe_SiO3%con,probe_SiO3%atom,success )
!print *,'##################### 1 ###############'
!         if(success)then
!            write(77,*) kmc_step,ib,rxyz(ia_OH,3)
!         end if
!      case(2) ! reaction 2
!         ia_OSiOH3 = irxn(2,ib,1)
!         booltmp = (ib < itop_nonp_SiO4)
!         call water_attack(booltmp,ia_OSiOH3,success)
!print *,'##################### 2 ###############'
!      case(3)
!!         ia_OH = irxn(3,ib,1)
!!         iOH  =  irxn(3,ib,2)
!         call condensation(ia_OH,iOH)
!print *,'##################### 3 ###############'
!      case(4)
!         ia_OSiO = irxn(4,ib,1)
!         booltmp = (ib < itop_nonp_SiO4)
!         call water_attack(booltmp,ia_OSiO,success)
!print *,'##################### 4 ###############'
!      end select
!!
!!     advance the kmc clock
!      delta_time = -log(rand())/rate_sum
!      time = time + delta_time
!!
!      call cpu_time(timep(2))
!      call new_vlist
!!
!      call RELAX(NRELAX, 1, natom, kBoltzT)
!
!      call cpu_time(timep(3))
!
!!     call relax_bonds( count(atom(1:natom) == iOxygen) )
!      call relax_bonds( min(10,count(atom(1:natom) == iOxygen)) )
!
!
!      call cpu_time(timep(4))
!      call date_and_time(date,ctime)
!!
!      write(* ,'(i7,a10,2x,a10,6(3x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
!         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
!      write(25,'(i7,a10,2x,a10,6(2x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
!         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
!      top_atom = maxval(rxyz(1:natom,3))
!      write(24,'(f12.6,2i7,f12.6)') time,kmc_step,natom,top_atom
!!
!      if (make_movie) then
!         call new_frame_file(imve,'frame',kmc_step)
!         call write_frame(imve,1,natom)
!      end if
!!
!      if ( mod(kmc_step,ipconfig) == 0 .or. (kmc_step == ntotal) ) then
!!--------Periodically check the connectivity and write out xyz coordinates,
!!--------RDF, Connectivity, etc.
!!         call check_proximity(ok,i)
!!         if (.NOT.ok) then
!!            write(*,*) 'proximity array is NOT consistent'
!!            write(*,*) 'check_proximity ',ok,i
!!            stop
!!         end if
!!         call check_proximity2(1,natom,ok,i)
!!         if (.NOT.ok) then
!!           write(*,*) 'proximity array is NOT consistent'
!!           write(*,*) 'check_proximity2 ',ok,i
!!           stop
!!         end if
!         call write_xyz(15)
!         rewind(15)
!         call write_proximity(20)
!         rewind(20)
!         call write_rand_state(21)
!         rewind(21)
!         call den_profile_z(1,natom,(5.0_wp/angstrom))
!         call write_den_profile_z(22,kmc_step)
!      end if
!
!      end do main_loop
!!----------------------------------------------------------------------
!!
!      write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
!      write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
!      write(*,*) 'N_OH= ',count(atom(1:natom) == iOxygenH)
!!
!      close(22)
!      OPEN(unit= 22,file= 'analysis_dens.out')
!      call den_profile(1,natom,(10.0_wp/angstrom))
!      call write_den_profile_z(22,kmc_step)
!      close(22)
!
!!      close(23)
!!      OPEN(unit= 23,file= 'analysis_void_SiO4.out')
!!      call voidage_profile(probe_SiO4,dlayer)
!!      call write_voidage(23)
!!      close(23)
!
!      zl =  5.0_wp/angstrom
!      zu = 25.0_wp/angstrom
!      vo2 = voidage_calc(50000,zl,zu,probe_O2)
!      vn2 = voidage_calc(50000,zl,zu,probe_N2)
!      vco2 = voidage_calc(50000,zl,zu,probe_CO2)
!      open(unit=40,file='analysis_void_5-25.out')
!      write(*,*) zl,zu
!      write(*,*) 'voidage O2 ',vO2
!      write(*,*) 'voidage N2 ',vN2
!      write(*,*) 'voidage CO2 ',vCO2
!      write(40,*) zl,zu
!      write(40,*) 'voidage O2 ',vO2
!      write(40,*) 'voidage N2 ',vN2
!      write(40,*) 'voidage CO2 ',vCO2
!      close(40)
!!
!      zl =  5.0_wp/angstrom
!      zu = 15.0_wp/angstrom
!      vo2 = voidage_calc(50000,zl,zu,probe_O2)
!      vn2 = voidage_calc(50000,zl,zu,probe_N2)
!      vco2 = voidage_calc(50000,zl,zu,probe_CO2)
!      open(unit=41,file='analysis_void_5-25.out')
!      write(*,*) zl,zu
!      write(*,*) 'voidage O2 ',vO2
!      write(*,*) 'voidage N2 ',vN2
!      write(*,*) 'voidage CO2 ',vCO2
!      write(41,*) zl,zu
!      write(41,*) 'voidage O2 ',vO2
!      write(41,*) 'voidage N2 ',vN2
!      write(41,*) 'voidage CO2 ',vCO2
!      close(41)
END PROGRAM porous_silica

