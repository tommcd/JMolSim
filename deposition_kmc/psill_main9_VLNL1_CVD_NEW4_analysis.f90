
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
INCLUDE 'netcdf_array_mod.f90'
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
!INCLUDE 'coordinates_mod.f90'
INCLUDE 'coordinates_mod_2DP.f90'
INCLUDE 'connectivity_mod.f90'
INCLUDE 'check_structure_mod_Silica.f90'
!INCLUDE 'check_structure_mod_TEOS.f90'
INCLUDE 'teos_atomlst_mod.f90'
INCLUDE 'nlist_mod3.f90'
INCLUDE 'verlet_list_nl.f90'
INCLUDE 'lattice_util_mod.f90'
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
!INCLUDE 'deposition5_mod.f90'
INCLUDE 'relax_VL_mod.f90'
INCLUDE 'select_atom_mod.f90'
!INCLUDE 'attach_tetrahedra_mod.f90'
INCLUDE 'attach_molecule_mod.f90'
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
INCLUDE 'kmc_reactions_mod_silica.f90'


PROGRAM porous_silica
      USE precision_mod
      USE command_line_mod
      USE command_arg_mod
      USE phys_const_mod
      USE global_vars_mod
      USE constants_mod
      USE atom_types_mod
      USE list_mod
      USE listset_mod
!     USE sets
      USE netcdf_arrays
      USE grid_data_mod
      USE connectivity_mod
      USE check_structure_mod
      USE coordinates_mod
      USE rand_mod
      USE frames_mod
      USE Relax_mc_mod
      USE nlist_mod
      USE verlet_list_mod
      USE HKNonLattice_mod
!     USE hk_lattice_mod
!     USE hk_lattice_nonlattice_mod
      USE unionfind_mod
      USE hk_lattice_uf_mod
      USE percolation_lattice_mod
      USE density_profile_mod
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
      USE kmc_reactions_mod
      USE lattice_util_mod
      implicit none
      logical:: make_movie
      real(wp):: top_atom,void_crit,rverlet,dlayer,Hmax
      real(wp):: dLflexible, dLfixed, dlkeep,quat(4),aa(3,3),ulj,uel,r1,rc,duel
      real(wp):: Rbox_min(3),Rbox_max(3),Rbox_top(3),del_nlist,fover_henry,del_bin
      real(wp):: z_top,dens_gcm3,volbox,rcutH,dl,Uljel_1,Khenry,facc,qj,sumqj,rv(3),energy
      integer:: i,itmp,j,ipconfig,ntotal,imve,narg,nc,ii,k,it,ix,iy,iz,stat,lb_rxn,movie_freq
      real:: timep(0:10)
      character:: date*8,ctime*10
      character(len=32):: ctmp,c6*6,c5*5,c32
      character(len=132):: coordfile,restart_coordfile,carg,ftyp*3,ctrlfile
      integer:: io_config, io_config_final, io_topol, io_ctrl, io_restart_xyz
      integer:: n_try_void,n_try_Henry
      logical:: success,ok,booltmp,restart
      real(wp),allocatable:: massi(:,:),mass(:,:)
      integer,allocatable:: listOH(:),listALL(:),listSph(:)
      integer:: nSph,noxdel,n_oh,n_si,n_ob
      type(probe_mol):: pr,pr1
      type(statistic):: s_lj,s_el,s_tot
      type(grid3D):: Grid
integer,parameter:: n_vapor_species = 2
integer,parameter:: iSiOH4 = 1
!integer,parameter:: iTeos = 1
integer,parameter:: iH2O = 2
      ! density, voidage, Henry's law constants
      real(wp),allocatable:: nden(:,:,:,:),ndenz(:,:)
      real(wp),allocatable:: voidage(:,:,:,:)
      real(wp),allocatable:: henryk(:,:,:,:)
      logical,allocatable:: Accessible(:,:,:,:),non_reacting(:,:,:)
      integer,allocatable:: Lattice(:,:,:,:)
integer:: ncluster(n_vapor_species)
type(mylist)::NumAcc(n_vapor_species)
integer:: nbinmax,nz_top_atom,nz_top,ib(3),irxn,ibin
integer:: nx,ny,nz,ndel,idel(1)
real(wp),allocatable:: rate(:,:,:,:),ndenrxn(:,:,:,:)
real(wp):: rate_sum, psum, volbin, HK_H2O,HK_SiOH4
type(mylist):: fixed_list
type(mylist):: atmlst(3),oh_lst,ohb_lst
integer:: ia_OH
integer:: ia_OSiOH3
integer:: ia_OSiO
integer:: ir,iOH
real(wp):: time,delta_time,xi1
real(wp),parameter:: rohbl = 4.0_wp/Angstrom
!real(wp),parameter:: volohbl = (4.0_wp/3.0_wp)*pi*(rohbl**3)
integer:: n_neigh,neighbors(2000)
real(wp):: P_SiOH4,P_H2O,T_Kelvin,zu,zl,vco2,vo2,vn2
logical:: H2O_present
logical:: hydrothermal = .false.
!
integer:: n_oh_pairs,iSi1,iSi2
real(wp):: rv_SiSi(3),r_SiSi
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
      read(io_ctrl,*) c32,void_crit               ; write(*,*) trim(c32), void_crit
      read(io_ctrl,*) c32,P_SiOH4                 ; write(*,*) trim(c32), P_SiOH4
      read(io_ctrl,*) c32, del_nlist              ; write(*,*) trim(c32), del_nlist
      read(io_ctrl,*) c32, rverlet                ; write(*,*) trim(c32), rverlet
      read(io_ctrl,*) c32, NRELAX                 ; write(*,*) trim(c32), NRELAX
      read(io_ctrl,*) c32, Hmax                   ; write(*,*) trim(c32), Hmax
      read(io_ctrl,*) c32, ipconfig               ; write(*,*) trim(c32), ipconfig
      read(io_ctrl,*) c32, make_movie,movie_freq  ; write(*,*) trim(c32), make_movie, movie_freq
      read(io_ctrl,*) c32, Hydrothermal           ; write(*,*) trim(c32), Hydrothermal
      read(io_ctrl,*) c32, n_try_void             ; write(*,*) trim(c32), n_try_void
      read(io_ctrl,*) c32, n_try_Henry            ; write(*,*) trim(c32), n_try_Henry
      read(io_ctrl,*) c32, fover_henry            ; write(*,*) trim(c32), fover_henry
      rcutH = 10.0_wp/Angstrom
      del_nlist = del_nlist/Angstrom
      rverlet = rverlet/Angstrom
      Hmax = Hmax/angstrom
      close(io_ctrl)
!
      kBoltzT = K_ev*T_Kelvin
      P_SiOH4 = P_SiOH4/upressure
      P_H2O = fwater*P_SiOH4
      if(hydrothermal) then
         P_H2O = P_SiOH4
         P_SiOH4 = 0.0_wp
      end if
      H2O_present = (P_H2O > tiny(0.0))
!
      print '(a,g18.9,a)','ulength   = ',ulength,' m'
      print '(a,g18.9,a)','uenergy   = ',uenergy,' J'
      print '(a,g20.10,a)','upressure = ',upressure,' Pa'
      print '(a,g18.9,a,g18.9,a)','P_SiOH4 = ',P_SiOH4,' = ',P_SiOH4*upressure,' Pa'
      print '(a,g18.9,a)',        '        = ',P_SiOH4*upressure/torr2Pa,' torr'
      print '(a,g18.9,a)','kT   = ',kBoltzT,' eV'
      print *,'T = ',T_Kelvin,' K'

!
! Read in the initial molecular system file header
      select case(ftyp)
      case('xyz')
         call read_xyz_header(io_config,natom,boxl)
         call get_free_file_unit(io_topol)
         open(unit=io_topol,file=trim(coordfile(1:nc - 3))//'con')
      case('pdb')
         call read_pdb_header(io_config,natom,boxl)
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select
!
      !natom_max = natom
      natom_max = natom + 7*ntotal  ! 3d
boxl(3) = max(boxl(3),2.0_wp*Hmax)
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
      allocate(charge(natom_max))
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
      write(*,'(a,f14.8)') 'qo = ',qi(iOxygen)
      write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)

      call count_atoms(1,natom,atom,n_atomtype)
      call write_atoms(6,n_atomtype)
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

!call Init_HKNonLattice(natom_max)
!call HKNonLattice(natom,proximity,n_cluster,atomL)
!n_cluster_old = n_cluster
!write(*,*) 'n_cluster = ',n_cluster
!do i = 1,natom
!   if (ANY(rxyz(i,:)>=boxl2).or.ANY(rxyz(i,:)<=-boxl2)) print '(A,I6,3F12.6)','BOUNDARY ', i,rxyz(i,:)
!end do
!Rbox_min = (/ -boxl2(1), -boxl2(2),-boxl2(3)-5.0_wp/Angstrom /)
!Rbox_max = (/  boxl2(1),  boxl2(2), boxl2(3)+5.0_wp/Angstrom /)
!call init_density_profile(Rbox_min,Rbox_max,5.0_wp/Angstrom,ntyp,nden,ndenz)
!call den_profile(1,natom,Rbox_min,Rbox_max,5.0_wp/Angstrom,nden,nndenz)
!call write_den_profile_z(6,0,Rbox_min(3),5.0_wp/Angstrom,ndenz)
!call den_profile_z(1,natom,Rbox_min(3),Rbox_max(3),5.0_wp/Angstrom,ndenz)
!call write_den_profile_z(6,0,Rbox_min(3),5.0_wp/Angstrom,ndenz)
!call print_mat(6,int(nden(:,:,:,1)))
!call print_mat(6,int(nden(:,:,:,2)))
!call print_mat(6,int(nden(:,:,:,3)))
!do i = 1,ntyp
!print *,i,atom_name(i),maxval(nden(:,:,:,i))
!end do
!allocate(Accessible(7,7,8,n_vapor_species))
!call INIT_NLIST(Rbox_min,Rbox_max,del_nlist,natom_max)
!if (del_nlist < maxval(sigLJ)) STOP 'Neighbor list is too small'
!call NEW_NLIST(1,natom)
!call  PRUNE_BLOCKED_CELLS(sum(nden,dim=4),0.035_wp,Accessible(:,:,:,1))
!print *,'count(Accessible(:,:,:,1)) = ',count(Accessible(:,:,:,1))
!print *,'count(sum(nden,dim=4) > 0.0_wp) = ',count(sum(nden,dim=4) > 0.0_wp)
!deallocate(nden,ndenz,Accessible)

! Initially the box should be centered at the origin.
      print '(a,3f12.6)','boxl old = ',boxl*Angstrom

! Set the size of cell bins, the depth of the substrate below z=0 plane to keep
! and the flexible and fixed regions of it
      del_bin = del_nlist   ! 5.0_wp/Angstrom
      dLflexible = -5.0_wp/Angstrom
      !dLfixed = 5.0_wp/Angstrom
      !dLkeep = dLflexible - dLfixed
      dLkeep = min(minval(rxyz(:,3)),-15.0_wp/angstrom)

      print *,' dLflexible = ', dLflexible*Angstrom
      !print *,'    dLfixed = ', dLfixed*Angstrom
      print *,'     dLkeep = ', dLkeep*Angstrom

!      if (kmc_step == 0) then
!         ! translate the box by -L_z/2
!         rxyz(1:natom,3) = rxyz(1:natom,3) - boxl2(3)
!   
!         ! and set the new Lz to Max(2*Hmax,2*Lz)
!         boxl(3) = max(2.0_wp*boxl(3),2.0_wp*Hmax)
!         call SET_BOX_SIZE(boxl)
!   
!         ! Delete all atoms below dLkeep
!         i = 0
!         do
!            i = i + 1
!            if(i > natom) exit
!            if (rxyz(i,3) < dlkeep) then
!               call delete_atom(i)   !; print *,'delete ',i
!               i = i - 1
!            end if
!         end do
!      end if

! Make a list of the fixed atoms
      call new_list(fixed_list,natom)
      do i = 1,natom
         if (rxyz(i,3) < dLflexible) call add_to_list(fixed_list, i)
      end do
      ! and renumber the fixed atoms to go from 1,2,3,...,nfixed
      nfixed = fixed_list%n
      do i = 1,nfixed
         call swap_atom(i,fixed_list%ind(i))   !; print *,'swap ',i,fixed_list%ind(i)
      end do
      if (nfixed /= count( rxyz(:,3) < dLflexible )) stop 'nfixed is wrong'
      call delete_list(fixed_list) ! finished with this list

!! Can also do it this way
!      ! Renumber the fixed atoms to go from 1,2,3,...,nfixed
!      nfixed = count( rxyz(:,3) < dLflexible )
!      do i = 1,natom
!         if (rxyz(i,3) < dLflexible) cycle
!         do j = i+1, natom
!            if (rxyz(j,3) < dLflexible) then
!               call swap_atoms(i,j)   !; print *,'swap ',i,j
!               exit
!            end if
!         end do
!      end do

      print *,'nfixed = ',nfixed
      print *,' natom = ',natom

! Set the upper & lower corners of the cuboid which defines the reaction region
! & regionwhere the neighbor list is defined.
      ! Rbox_min = -boxl2; Rbox_min(3) = Rbox_min(3)-5.0_wp/Angstrom
     ! Rbox_min = (/ -boxl2(1), -boxl2(2), dLflexible /)
      Rbox_min = (/ -boxl2(1), -boxl2(2), dLKeep /)
      Rbox_max = (/  boxl2(1),  boxl2(2), boxl2(3) /)

print '(a,3f12.6)','Rbox_min = ',Rbox_min*Angstrom
print '(a,3f12.6)','Rbox_max = ',Rbox_max*Angstrom


! !atom(1:nfixed) = iOw
! call new_frame_file(imve,'frame',20000,'xyz')
! call write_xyz_file(imve,1,natom,atom,rxyz,boxl, &
!                     box_corners=cuboid( Rbox_min,Rbox_max ), &
!                     comment=8//' Box corners added')
! call new_frame_file(imve,'frame',20000,'con')
! call write_conn(imve,1,natom,proximity,'connectivity at step = '//0)




! Define the 3D Grid    
      call new_grid3d(Rbox_min,Rbox_max,del_bin,Grid)

      print '(a,3f12.6)','Grid cell sides = ',Grid%del
      print '(a,3i7)',   'Grid cell bins  = ',Grid%nbin
      Nx = Grid%nbin(1)
      Ny = Grid%nbin(2)
      Nz = Grid%nbin(3)
      nbinmax = Nx*Ny*Nz
      volbin = Grid%volbin

      allocate(nden(nx,ny,nz,ntyp),ndenz(Nz,ntyp))
      allocate(voidage(nx,ny,nz,n_vapor_species))
      allocate(HenryK(nx,ny,nz,n_vapor_species))
      allocate(Accessible(nx,ny,nz,n_vapor_species))
      accessible = .FALSE.
      allocate(Lattice(nx,ny,nz,n_vapor_species),non_reacting(nx,ny,nz))
      do i = 1,n_vapor_species
         call new_list(NumAcc(i),Nz)
      end do
      allocate(rate(nrxn,nx,ny,nz))
      allocate(rxns(nx,ny,nz))
      allocate(ndenrxn(nx,ny,nz,natom_type_rxn))

! Initialize and create the neighbour list
      call INIT_NLIST(Rbox_min,Rbox_max,del_nlist,natom_max)
      if (del_nlist < maxval(sigLJ)) STOP 'Neighbor list is too small'
      call NEW_NLIST(1,natom)

! Initialize and create the Verlet list
      call Init_Verlet_List(rv0=rverlet,rcut0=2.6_wp/Angstrom,nat_max=natom_max)
      call new_vlist(nfixed+1,natom)

!
! Initialize both the lattice and non-lattice versions of the Hoshen-Kopelman
! cluster labeling algorithms and allocate the percolation analysis arrays.
!
      call Init_HKNonLattice(natom_max)
      call HKNonLattice(natom,proximity,n_cluster,atomL)
      n_cluster_old = n_cluster
      write(*,*) 'n_cluster = ',n_cluster
!      call hk_initialize((nx*ny*nz)/2)


!
! output files
!      OPEN(unit= 15,file= 'psil_config.pdb')
!      OPEN(unit= 20,file= 'psil_topology.out')
!      OPEN(unit= 21,file= 'psil_rand.out')
!      OPEN(unit= 22,file= 'psil_den_z.out')
!      OPEN(unit= 23,file= 'psil_void_SiO4.out')
!      OPEN(unit= 24,file= 'psil_topatom.out')
!      OPEN(unit= 25,file= 'psil_time.out')
!      OPEN(unit= 29,file= 'psil_kmc_history.out')
!      OPEN(unit= 77,file= 'reaction_01.out')
!      OPEN(unit= 78,file= 'reaction_02.out')
!      OPEN(unit= 80,file= 'reaction_03.out')
!      OPEN(unit= 81,file= 'reaction_04.out')
!      OPEN(unit= 82,file= 'psil_bond_switch.out')
!
!      do i = 1,size(atmlst)
!         call new_list(atmlst(i),1750)
!      end do
!      call new_list(oh_lst,100)
!      call new_list(ohb_lst,100)
!      call cpu_time(timep(0))
!
!----------------------------------------------------------------------
!      time = 0.0_wp
!      main_loop: do while (kmc_step < ntotal)
!      call cpu_time(timep(1))
!      kmc_step = kmc_step + 1

      call assign_charge(1,natom)

      top_atom = maxval(rxyz(1:natom,3))
      nz_top_atom = int((top_atom-Grid%rlower(3))*Grid%deli(3)) + 1
      nz_top = nz_top_atom + 2
      z_top = nz_top*Grid%del(3)
      Rbox_top = (/ boxl2(1), boxl2(2), z_top /)
print *,'z_top = ',z_top*Angstrom

      call den_profile(nfixed+1,natom,Grid,nden,ndenz)

!! Calculate the voidage in the cells up to z-layer=nz_top
!     !call voidage_profile(probe_teos, Rbox_min, Rbox_max, del_bin, n_try_void, voidage)
!      call voidage_profile( probe_SiO4, Rbox_min, Rbox_top, del_bin, n_try_void, voidage(:,:,:,iSiOH4) )
!      if (H2O_present) call voidage_profile( probe_tip3p,Rbox_min, Rbox_top, del_bin, n_try_void, voidage(:,:,:,iH2O) )
!
!! Mark the cells with voidage > critical value for each Vapor phase reactant
!! (e.g. for both Si(OH)4 or TEOS and H2O)
!
!      where (voidage(:,:,1:nz_top,:) >= void_crit)
!             Lattice(:,:,1:nz_top,:) = 1
!      elsewhere
!             Lattice(:,:,1:nz_top,:) = 0
!      end where

      call voidage_check_lattice( probe_SiO4, Rbox_min, Rbox_top, del_bin, void_crit, n_try_void, Lattice(:,:,:,iSiOH4) )
      call voidage_check_lattice( probe_tip3p,Rbox_min, Rbox_top, del_bin, void_crit, n_try_void, Lattice(:,:,:,iH2O) )

!
! Percolation analysis
      do i = 1, n_vapor_species
         ! Use the Hoshen-Kopelman algorithm to label the clusters of connected cells for
         ! each vapor phase reactant
         call hk_lattice_uf( Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i) )
         ! Alternative versions of algorithm (for comparison)
         ! call hk_lattice_3d(Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i))
         ! call hk_lattice_3dnl(Lattice(:,:,1:Nz_top,i),box_periodic,ncluster(i))
         !
         ! Next mark the cell which are accessible (from the vapor phase)
         call ACCESSIBLE_REGIONS_Z( Lattice(:,:,1:Nz_top,i), 1, Nz_top, &
                                 Accessible(:,:,1:Nz_top,i), NumAcc(i) )
print *,i,ncluster(i)
      end do


      do j = 1, n_vapor_species
         print *,'# NumAcc%n = ',NumAcc(j)%n
         do i = 1,NumAcc(j)%n
            print *,i,NumAcc(j)%ind(i)
         end do
      end do
call cpu_time(timep(2))
print *,'time for voidage check = ',timep(2)-timep(1)

! calculate the Henry's law constant only in the accessible cells
    ! call Henrys_law_profile(probe,rlower,rupper,ntrial,fover_henry,rcutH,dl,nb,grid)
      call Henrys_law_grid(probe_SiO4,Accessible(:,:,1:Nz_top,iSiOH4),grid,Nz_top, &
                           n_try_henry,fover_henry,rcutH,henryk(:,:,1:Nz_top,iSiOH4))
!      call Henrys_law_grid(probe_TEOS,Accessible(:,:,1:Nz_top,iTEOS),grid,Nz_top, &
!                           n_try_henry,fover_henry,rcutH,henryk(:,:,1:Nz_top,iTEOS)
      if (H2O_present) call Henrys_law_grid(probe_tip3p,Accessible(:,:,1:Nz_top,iH2O),grid,Nz_top, &
                           n_try_henry,fover_henry,rcutH,henryk(:,:,1:Nz_top,iH2O))

call cpu_time(timep(3))
print *,'time for Henrys_law_grid H2O= ',timep(3)-timep(2)

      non_reacting = .false.
      lb_rxn = 3  !abs(dLflexible)/Grid%del(3)
!print *, 'abs(dLflexible) = ',abs(dLflexible)
!print *, 'Grid%del(3) = ',Grid%del(3)
!print *, 'lb_rxn =  ',lb_rxn
      non_reacting(:,:,1:lb_rxn) = .true.

      call den_profile_rxn(nfixed+1,natom,Grid,ndenrxn)


!      rate_sum = 0.0_wp
!      rate = 0.0_wp
!      forall(i=1:4)rxns%r(i)%n=0
!      forall(i=1:4,j=1:2)rxns%r(i)%atom(j)=0
!
!      do iz = lb_rxn, Nz_top
!      do iy = 1, Ny
!      do ix = 1, Nx
!         ibin = icell(ix,iy,iz)
!         !print *,'ibin = ',ibin,icell(ix,iy,iz)
!
!         call NEIGCELL(ibin,1,n_neigh,neighbors)
!
!         call list_atoms_in_cells(ibin,atmlst)
!
!!         call list_atoms_in_cells(n_neigh,neighbors,atmlst)
!
!         call rand_from_list(atmlst(L_OH),ia_OH)
!         call rand_from_list(atmlst(L_OSiOH3),ia_OSiOH3)
!         call rand_from_list(atmlst(L_OSiO),ia_OSiO)
!
!k_1f = A_1f*exp(-E_A1f/kBoltzT)
!k_1r = A_1r*exp(-E_A1r/kBoltzT)
!k_2r = A_2r*exp(-E_A2r/kBoltzT)
!
!
!         if (ia_OH /= 0) then
!!
!!                         k_1F
!! (1)  |-OH  +  Si(OH)4   --->   |-Si(OH)3  +  H2O
!!
!            k_1f = A_1f*exp(-E_A1f/kBoltzT)
!            !HK_SiOH4 = HenryK(ix,iy,iz,iSiOH4)
!            HK_SiOH4 = HenryK_Local_Average(henryk(:,:,:,iSiOH4),ix,iy,iz)
!!            print *,'HenryK(ix,iy,iz,iSiOH4) = ',HenryK(ix,iy,iz,iSiOH4)
!!            print *,'               HK_SiOH4 = ',HK_SiOH4
!            rate(1,ix,iy,iz) = HK_SiOH4*P_SiOH4*ndenrxn(ix,iy,iz,iOxygenH)
!            if (.NOT.Accessible(ix,iy,iz,iSiOH4)) rate(1,ix,iy,iz) = 0.0_wp
!            rate_sum = rate_sum + rate(1,ix,iy,iz)
!            rxns(ix,iy,iz)%r(1)%n = 1
!            rxns(ix,iy,iz)%r(1)%atom(1) = ia_OH
!         end if
!
!
!         if (H2O_present .and. (ia_OSiOH3 /= 0)) then
!!
!!                         k_1R
!! (2)  |-Si(OH)3  +  H2O  --->   |-OH  +  Si(OH)4
!!
!            !HK_H2O = henryk(ix,iy,iz,iH2O)
!            HK_H2O = HenryK_Local_Average(henryk(:,:,:,iH2O),ix,iy,iz)
!            rate(2,ix,iy,iz) = (k_1r/k_1f)*HK_H2O*P_H2O*ndenrxn(ix,iy,iz,iOSiOH3)
!            if (.NOT.Accessible(ix,iy,iz,iSiOH4)) rate(2,ix,iy,iz) = 0.0_wp
!            rate_sum = rate_sum + rate(2,ix,iy,iz)
!            rxns(ix,iy,iz)%r(2)%n = 1
!            rxns(ix,iy,iz)%r(2)%atom(1) = ia_OSiOH3
!         end if
!
!
!         if (ia_OH /= 0) then
!!
!!                         k_2F
!! (3)     |-OH  +  HO-|   --->   |-O-|  +  H2O
!!
!            call find_atoms(ia_OH, nfixed+1, natom, iOxygenH, rohbl, oh_lst%n, oh_lst%ind)
!            ohb_lst%n = 0
!            ohb_lst%ind = 0
!            do i = 1,oh_lst%n
!               iOH = oh_lst%ind(i)
!!print *
!!print *,'ia_OH = ',ia_OH
!!call write_atom_info(6,ia_OH)
!!print *,'iOH = ',iOH
!!call write_atom_info(6,iOH)
!               if (OH_groups_OK(ia_OH,iOH)) call add_to_list(ohb_lst,iOH)
!            end do
!            n_oh_pairs = ohb_lst%n
!            call rand_from_list(ohb_lst,iOH)
!            if (iOH /= 0) then
!               call find_1st_atomtype(ia_OH,iSilicon,iSi1)
!               call find_1st_atomtype(iOH,iSilicon,iSi2)
!               rv_SiSi = rxyz(iSi1,:) - rxyz(iSi2,:)
!               call pbc(rv_SiSi)
!               r_SiSi = sqrt(dot_product(rv_SiSi,rv_SiSi))
!               k_2f = A_2f*exp(-E_A_2F(r_SiSi*Angstrom)/kBoltzT)
!               !
!               ! v_3 = kT (k_3/k_1) N_{OH-OH}
!               !
!               rate(3,ix,iy,iz) = kBoltzT*(k_2f/k_1f)*n_oh_pairs
!               if (.NOT.Accessible(ix,iy,iz,iH2O)) rate(3,ix,iy,iz) = 0.0_wp
!               rate_sum = rate_sum + rate(3,ix,iy,iz)
!               rxns(ix,iy,iz)%r(3)%n = 2
!               rxns(ix,iy,iz)%r(3)%atom(1) = ia_OH
!               rxns(ix,iy,iz)%r(3)%atom(2) = iOH
!            end if
!         end if
!
!
!         if (H2O_present .and. (ia_OSiO /= 0)) then
!!
!!                         k_2R
!! (4)     |-O-|  +  H2O   --->   |-OH  +  HO-|
!!
!print *,'lb_rxn = ',lb_rxn
!print *,'Nz_top = ',Nz_top
!print *,ix,iy,iz
!            !HK_H2O = henryk(ix,iy,iz,iH2O)
!            HK_H2O = HenryK_Local_Average(henryk(:,:,:,iH2O),ix,iy,iz)
!            rate(4,ix,iy,iz) = (k_2r/k_1f)*HK_H2O*P_H2O*ndenrxn(ix,iy,iz,iOSiO)
!            if (.NOT.Accessible(ix,iy,iz,iH2O)) rate(4,ix,iy,iz) = 0.0_wp
!            rate_sum = rate_sum + rate(4,ix,iy,iz)
!            rxns(ix,iy,iz)%r(4)%n = 1
!            rxns(ix,iy,iz)%r(4)%atom(1) = ia_OSiO
!         end if
!
!
!      end do
!      end do
!      end do
!!
!      xi1 = rand()
!      psum = 0.0_wp
!      bin_loop: do iz = lb_rxn, Nz_top
!      do iy = 1, Ny
!      do ix = 1, Nx
!      do ir = 1, nrxn
!         psum = psum + rate(ir,ix,iy,iz)/rate_sum
!         if (psum >= xi1) then
!            ib = (/ ix,iy,iz /)
!            irxn = ir
!            exit bin_loop
!         end if
!      end do
!      end do
!      end do
!      end do bin_loop
!
!      write(*, '(g18.8,i4,3i9)') rate_sum,irxn,ib
!      write(29,'(g18.8,i4,3i9)') rate_sum,irxn,ib
!!
!!
!print *,'irxn = ',irxn
!print *,'ib   = ',ib
!
!!
!!     carry out reaction irxn in bin ib
!      select case(irxn)
!      case(1) ! reaction 1
!         ia_OH = rxns(ib(1),ib(2),ib(3))%r(1)%atom(1)
!         ndel = 0
!         idel(1) = 0
!         !call attach_tetrahedra(ia_OH,aSiOSi,aSiO,overlap_check=.false.,success=success)
!         !call attach_mol( ia_OH,ndel,idel,bondl_SiO,angle_SiOSi,4,transpose(probe_SiO3%r), &
!         !                    probe_SiO3%conect,probe_SiO3%atom,success )
!         call attach_mol( ia_OH,ndel,idel,bondl_SiO,angle_SiOSi,4,transpose(probe_SiO3%r), &
!                             probe_SiO3%conect,probe_SiO3%atom,force_attachment=.true. )
!         success=.true. ! because we forced it
!         if (success) then
!            print *,'attach in cell ib = ',ib
!print *,'rxns(ib(1),ib(2),ib(3))%r(1)%atom',rxns(ib(1),ib(2),ib(3))%r(1)%atom(:)
!print *,'ia_OH = ',ia_OH
!            write(*,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OH,3)*Angstrom
!            write(77,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OH,3)*Angstrom
!if(rxyz(ia_OH,3)*Angstrom<-5.0_wp)stop 'rxyz(ia_OH,3)*Angstrom<-5.0_wp'
!         end if
!      case(2) ! reaction 2
!         ia_OSiOH3 = rxns(ib(1),ib(2),ib(3))%r(2)%atom(1)
!         booltmp = .NOT.Accessible(ib(1),ib(2),ib(3),iSiOH4)
!         call water_attack(booltmp,ia_OSiOH3,success)
!         if (success) then
!            write(*,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OSiOH3,3)*Angstrom
!            write(78,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OSiOH3,3)*Angstrom
!         end if
!      case(3)
!         ia_OH = rxns(ib(1),ib(2),ib(3))%r(3)%atom(1)
!         iOH  =  rxns(ib(1),ib(2),ib(3))%r(3)%atom(2)
!print *, ia_OH,iOH
!call write_atom_info(6,ia_OH)
!call write_atom_info(6,iOH)
!         call condensation(ia_OH,iOH)
!         write(*,*) 'H2O removal ',rxyz(ia_OH,3)*Angstrom
!         write(80,*) kmc_step,rxyz(ia_OH,3)*Angstrom
!      case(4)
!         ia_OSiO = rxns(ib(1),ib(2),ib(3))%r(4)%atom(1)
!         booltmp = .NOT.Accessible(ib(1),ib(2),ib(3),iSiOH4)
!         call water_attack(booltmp,ia_OSiO,success)
!         if (success) then
!            write(*,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OSiO,3)*Angstrom
!            write(81,'(i7,3i5,2x,f12.6)') kmc_step,ib,rxyz(ia_OSiO,3)*Angstrom
!         end if
!      end select
!
!!call check_proximity(1,natom,ok,i)
!!if (.NOT.ok) then
!!   write(*,*) 'proximity array is NOT consistent'
!!   write(*,*) 'check_proximity ',ok,i
!!   call write_atom_info(6,i)
!!   stop
!!end if
!!call check_proximity2(1,natom,atom,proximity,ok,i)
!!if (.NOT.ok) then
!!  write(*,*) 'proximity array is NOT consistent'
!!  write(*,*) 'check_proximity2 ',ok,i
!!  call write_atom_info(6,i)
!!  stop
!!end if
!
!
!
!
!!
!!     Advance the KMC clock
!      delta_time = -log(rand())/rate_sum
!      time = time + delta_time
!!
!      call cpu_time(timep(2))
!      call new_vlist(1,natom)
!!
!if ( mod(kmc_step,ipconfig) == 0 .or. (kmc_step == ntotal) ) then
!      call RELAX_MC(NRELAX, nfixed+1, natom, kBoltzT)
!end if
!
!      call cpu_time(timep(3))
!!
!!!     call relax_bonds( count(atom(1:natom) == iOxygen) )
!!      call relax_bonds( min(10,count(atom(1:natom) == iOxygen)) )
!
!
!      call cpu_time(timep(4))
!      call date_and_time(date,ctime)
!!
!      write(* ,'(i7,a10,2x,a10,6(3x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
!         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
!      write(25,'(i7,a10,2x,a10,6(2x,f0.3))') natom,date,ctime,timep(2)-timep(1), &
!         timep(3)-timep(2),timep(4)-timep(3),timep(4)-timep(0)
!!
!      top_atom = maxval(rxyz(1:natom,3))
!      write(24,'(f0.6,2i7,f14.6)') time,kmc_step,natom,top_atom
!!
!      if (make_movie.and.(mod(kmc_step,movie_freq)==0)) then
!         call new_frame_file(imve,'frame',kmc_step,'pdb')
!         call write_frame(imve,1,natom)
!      end if
!!
!      if ( mod(kmc_step,ipconfig) == 0 .or. (kmc_step == ntotal) ) then
!!--------Periodically check the connectivity and write out xyz coordinates,
!         call check_proximity(1,natom,ok,i)
!         if (.NOT.ok) then
!            write(*,*) 'proximity array is NOT consistent'
!            write(*,*) 'check_proximity ',ok,i
!            call write_atom_info(6,i)
!            stop
!         end if
!         call check_proximity2(1,natom,atom,proximity,ok,i)
!         if (.NOT.ok) then
!           write(*,*) 'proximity array is NOT consistent'
!           write(*,*) 'check_proximity2 ',ok,i
!           call write_atom_info(6,i)
!           stop
!         end if
!         call write_calc_den_profile_z(22,kmc_step,1,natom,5.0_wp/Angstrom)
!!--------RDF, Connectivity, etc.
!!         call write_xyz(15)
!!         rewind(15)
!!         call write_proximity(20)
!!         rewind(20)
!         call write_rand_state(21)
!         rewind(21)
!      end if
!!
!      end do main_loop
!!----------------------------------------------------------------------
!!
!      call count_atoms(1,natom,atom,n_atomtype)
!      call write_atoms(6,n_atomtype)
!      write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
!      write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
!      write(*,*) 'N_OH= ',count(atom(1:natom) == iOxygenH)
!!
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
CONTAINS


   PURE FUNCTION HenryK_Local_Average(henryk,ix,iy,iz)
      USE nlist_mod
      real(wp):: HenryK_Local_Average
      real(wp),intent(in):: henryk(:,:,:)
      integer,intent(in):: ix,iy,iz
      integer:: ic,inc,k,inx,iny,inz
      integer:: nn,neighbors(27),ncell
      real(wp):: hsum
      hsum = 0.0_wp
      ncell = 0
      ic = icell(ix,iy,iz)
      call neigcell(ic,1,nn,neighbors)
      do k = 1,nn
         inc = neighbors(k)
         !if(inc == ic) cycle
         call cell2xyz(inc,inx,iny,inz)
         ncell = ncell + 1
         hsum = hsum + henryk(inx,iny,inz)
      end do
      HenryK_Local_Average = hsum/ncell
   END FUNCTION HenryK_Local_Average

END PROGRAM porous_silica

