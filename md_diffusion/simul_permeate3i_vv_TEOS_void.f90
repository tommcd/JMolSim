
INCLUDE 'precision_mod.f90'
INCLUDE 'math_const_mod.f90'
INCLUDE 'phys_const_mod.f90'
INCLUDE 'command_line_mod.f90'
INCLUDE 'command_arg_mod.f90'
INCLUDE 'files_mod.f90'
INCLUDE 'rand_mod.f90'
INCLUDE 'sort_mod.f90'
INCLUDE 'grid_data_mod.f90'
!INCLUDE 'netcdf_array_mod.f90'
INCLUDE 'global_vars_mod_MD.f90'
INCLUDE 'HKNonLattice2.f90'
INCLUDE 'hk_lattice_mod.f90'
INCLUDE 'Keating_parameters_mod_vonAlfthan.f90'
INCLUDE 'Dreiding_parameters_mod.f90'
INCLUDE 'charmm_mod.f90'
INCLUDE 'seaton_mod.f90'
INCLUDE 'atom_types_mod_teos.f90'
INCLUDE 'constants_mod_teos.f90'
INCLUDE 'coordinates_mod.f90'
INCLUDE 'connectivity_mod.f90'
INCLUDE 'check_structure_mod_Teos.f90'
INCLUDE 'list_mod.f90'
INCLUDE 'teos_atomlst_mod.f90'
INCLUDE 'nlist_mod2.f90'
INCLUDE 'verlet_list_nl.f90'
INCLUDE 'find_atoms_mod.f90'
INCLUDE 'charges2_mod_teos.f90'
INCLUDE 'rotate_axis_mod.f90'
INCLUDE 'tetra_coords_mod.f90'
INCLUDE 'probe_mol_mod.f90'
INCLUDE 'bond_angle_types_mod_TEOS.f90'
INCLUDE 'bond_angle_list_mod.f90'
INCLUDE 'box_frames_mod.f90'
INCLUDE 'comvel_mod.f90'
INCLUDE 'vverlet_mod.f90'
INCLUDE 'glassbond_mod.f90'
INCLUDE 'repul_force_NL_vonAlfthan.f90'
INCLUDE 'quat2mat_mod.f90'
INCLUDE 'ran_point_sphere_mod.f90'
INCLUDE 'insert_mod_2.f90'
INCLUDE 'lj_el_mod_O2_N2_CO2.f90'
INCLUDE 'force_en_el_mod.f90'
INCLUDE 'tcf_mod.f90'
INCLUDE 'msd_mod3.f90'
INCLUDE 'vcf_mod3.f90'
INCLUDE 'rdf_mod.f90'
INCLUDE 'utility_mod.f90'
INCLUDE 'energy_keating_mod.f90'
INCLUDE 'force_keating_mod.f90'
INCLUDE 'statistics_mod.f90'
INCLUDE 'timer_mod.f90'
INCLUDE 'force_energy_minimize_mod.f90'
INCLUDE 'Steepest_Descent_mod.f90'
INCLUDE 'relax_mc_mod.f90'
INCLUDE 'attach_molecule_mod.f90'
INCLUDE 'voidage5_mod.f90'
INCLUDE 'Henrys_law_calc_mod.f90'


PROGRAM SIMBOX
      USE precision_mod
      USE phys_const_mod
      USE commandline_mod
      USE command_arg_mod
      USE precision_mod
      USE utility_mod
      USE files_mod
      USE constants_mod
      USE global_vars_mod
      USE coordinates_mod
      USE connectivity_mod
      USE list_mod
      USE bond_angle_list_mod
      USE rand_mod
      USE atom_types_mod
      USE Lj_el_mod
      USE box_frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE nlist_mod
      USE verlet_list_mod
      USE repul_energy_mod
      USE charges_mod
      USE force_en_el_mod
      USE insert_mod
      USE tcf_mod
      USE msd_mod
      USE vcf_mod
      USE rdf_mod
      USE Steepest_Descent_mod
      USE relax_mc_mod
      USE energy_keating_mod
      USE statistics_mod
      USE timer_mod
      USE find_atoms_mod
      USE probe_mol_mod
      USE check_structure_mod
      USE voidage_mod
      USE Henrys_law_calc_mod
      USE hk_lattice_mod
      implicit none
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      integer:: number_moves, nequib, nsamp, nprint_all, nfreq_movie, n_short
      logical:: make_movie
      integer:: ntau_msd, ntau_vcf
      real(wp):: del_nlist, fover,rl(NDIM),ru(NDIM),energy
      real(wp):: dt,qj,sumqj,T_kelvin,T_kinetic,sumx,sumy,sumz
      real(wp),allocatable:: mol_vec(:,:),mol_vec_h2o(:,:)
      real(wp):: EK,Ek_gas,Ek_O2,Ek_N2,Ek_CO2,Ek_tot,Ek_h2o
      real(wp):: dtreal,time,rr(NDIM),mtmp,ULJAtt,Ukeating,Ubond,Uang,uel
      real(wp),allocatable:: ULJEL(:)
      real(wp),allocatable::e_bond_CO2(:),e_angle_CO2(:),e_bond_gas(:),e_angle_gas(:)
      real(wp):: e_bond_teos,e_angle_teos
      real(wp):: e_bond_O2,e_bond_N2,Uold,Unew
      real(wp),allocatable:: massi(:,:),mass(:,:)
      real(wp),allocatable:: h2o_massi(:,:,:),h2o_mass(:,:,:),rad_h2o(:)
      real(wp),allocatable:: GAS_massi(:,:,:),GAS_mass(:,:,:),rad_GAS(:)
      real(wp),allocatable:: Ox_massi(:,:,:),Ox_mass(:,:,:)
      real(wp),allocatable:: N2_massi(:,:,:),N2_mass(:,:,:)
      real(wp),allocatable:: CO2_massi(:,:,:),CO2_mass(:,:,:)
      real(wp):: rad_Ox(3),rad_N2(3),rad_CO2(3)
      integer,allocatable:: mol_atoms(:), mol_atoms_h2o(:)
      integer,allocatable:: listOH(:),listALL(:),listSph(:),alst(:)
      integer:: nSph,nlst,nfixed,noxdel,n_oh,n_si,n_ob
      integer:: j,i,istep,stat,nc,ix,iy,iz,k
      integer:: itmp,narg,Ndf
      real(wp),allocatable:: coprxyz(:,:)
      real:: tp(0:10),sumtp(1:10) = 0.0
      character(len=32):: ctmp,c6*6,c5*5,c32
      character(len=132):: coordfile,restart_coordfile,restart_velfile,carg,ftyp*3,ctrlfile
      integer:: io_temp, io_energy, io_config, io_config_final, io_topol, io_ctrl
      integer:: io_restart_xyz, io_restart_v
      logical:: restart,ok,check_self_overlap=.true.
      type(tcf_type):: msd_gas, vcf_gas, msd_o2, vcf_o2, msd_n2, vcf_n2,msd_co2, vcf_co2, msd_h2o, vcf_h2o
      type(bond_angle_list) teos_bonds
      type(timer):: t(30)
      real(wp):: dLflexible, dLfixed, dlkeep,voidage(5),volbox,dens_gcm3,fdens_wt,rcutH,dl,Uljel_1,Khenry,facc
      integer:: nb(3),nlab,ntrial,nclusters
      integer,allocatable:: lab(:),Lattice(:,:,:)
      logical,allocatable:: Lb(:,:,:)
      real(wp),allocatable:: voidxyz(:,:,:),HenryK(:,:,:)
!
      print *,'Molecular Dynamics'
      print *,'Number of dimensions = ',NDIM
      print *,'Number of periodic dimensions = ',NPER
      if (NDIM < NPER) then
         print *,'The number of dimensions must be >= the number of periodic dimensions'
         STOP 'Change NDIM &/ NPER in coordinates_mod and recompile'
      end if
      call cpu_time(tp(0))
      call get_free_file_unit(io_energy)
      open(unit=io_energy,file='energy.out')
      call get_free_file_unit(io_temp)
      open(unit=io_temp,file='temperature.out')
!
      narg = command_argument_count()
      call next_command_argument(carg, "command",stat)
      if (narg < 3) then
         write(*,*)'usage :'
         write(*,*) trim(carg),'  control_file  atom_file  T[K]  {xyz_restart  vel_restart} '
         stop
      end if
!
      call next_command_argument(ctrlfile, "control file",stat)
      call get_free_file_unit(io_ctrl)
      open(unit=io_ctrl,file=trim(ctrlfile))
      read(io_ctrl,*) c32, number_moves        ; write(*,*) trim(c32), number_moves
      read(io_ctrl,*) c32, n_short             ; write(*,*) trim(c32), n_short
      read(io_ctrl,*) c32, nequib              ; write(*,*) trim(c32), nequib
      read(io_ctrl,*) c32, dtreal              ; write(*,*) trim(c32), dtreal
      read(io_ctrl,*) c32, nsamp               ; write(*,*) trim(c32), nsamp
      read(io_ctrl,*) c32, nprint_all          ; write(*,*) trim(c32), nprint_all
      read(io_ctrl,*) c32, make_movie          ; write(*,*) trim(c32), make_movie
      read(io_ctrl,*) c32, nfreq_movie         ; write(*,*) trim(c32), nfreq_movie
      read(io_ctrl,*) c32, ntau_msd            ; write(*,*) trim(c32), ntau_msd
      read(io_ctrl,*) c32, ntau_vcf            ; write(*,*) trim(c32), ntau_vcf
      read(io_ctrl,*) c32, n_O2,n_atom_o2      ; write(*,*) trim(c32), n_O2,n_atom_o2
      read(io_ctrl,*) c32, n_N2,n_atom_n2      ; write(*,*) trim(c32), n_N2,n_atom_n2
      read(io_ctrl,*) c32, n_CO2,n_atom_CO2    ; write(*,*) trim(c32), n_CO2,n_atom_CO2
      read(io_ctrl,*) c32, n_GAS, n_atom_gas   ; write(*,*) trim(c32), n_GAS, n_atom_gas
      read(io_ctrl,*) c32, n_h2o, n_atom_h2o   ; write(*,*) trim(c32), n_h2o, n_atom_h2o
      read(io_ctrl,*) c32, del_nlist           ; write(*,*) trim(c32), del_nlist
      read(io_ctrl,*) c32, fover               ; write(*,*) trim(c32), fover
      del_nlist = del_nlist/Angstrom

      call next_command_argument(coordfile, "input file",stat)
      nc = len_trim(coordfile)
print *,"'",trim(coordfile),"'"
print *,nc
      ftyp = coordfile(nc - 2:nc)
      call get_free_file_unit(io_config)
      open(unit=io_config,file=trim(coordfile))
      call next_command_argument(T_kelvin, "Temperature [K]",stat)
      kBoltzT = K_ev*T_kelvin
print *,'kBoltzT = ',kBoltzT, ' eV'
      call next_command_argument(restart_coordfile, "restart coord file",stat)
      if (stat /= 0) then
         restart = .false.
      else
         restart = .true.
         write (*,*) 'restarting from ', trim(restart_coordfile)
         call get_free_file_unit(io_restart_xyz)
         open(unit=io_restart_xyz,file=trim(restart_coordfile))
         call next_command_argument(restart_velfile, "restart vel file", stat)
         if (stat /= 0) then
            print *,'restart velocity file missing'
            stop
         else
            call get_free_file_unit(io_restart_v)
            open(unit=io_restart_v,file=trim(restart_velfile))
         end if
      end if
!
      dt = dtreal*dt_fs

      call init_probe_mols()

      rad_Ox = (/ sig_O_O2,sig_O_O2,0.0_wp /)*0.5_wp
      rad_N2 = (/ sig_n_N2,sig_n_N2,0.0_wp /)*0.5_wp
      rad_CO2 = (/ sig_C_CO2,sig_O_CO2,sig_O_CO2 /)*0.5_wp

!if(n_atom_gas/=33)stop 'n_atom_gas/=33'
n_atom_gas = 33
      allocate(mol_atoms(n_atom_gas))
      allocate(GAS_atom(1:n_GAS,n_atom_gas))
      allocate(GAS_charge(n_GAS,n_atom_gas))
      allocate(GAS_mass(1:n_GAS,n_atom_gas,3),GAS_massi(1:n_GAS,n_atom_gas,3))
      allocate(rad_GAS(n_atom_gas))
      !if(n_atom_gas==5 )mol_atoms = (/ iSilicon,iOxygenH,iOxygenH,iOxygenH,iOxygenH /)
      mol_atoms = pr_TEOS%atom
call new_bond_angle_list(TEOS_bonds,pr_TEOS%n)
call set_bond_angle_lists2(1,pr_TEOS%n,pr_TEOS%conect,pr_TEOS%atom,&
                           TEOS_bonds%nbond,TEOS_bonds%ibond,TEOS_bonds%kbond,TEOS_bonds%abond,&
                           TEOS_bonds%nang,TEOS_bonds%iang,TEOS_bonds%kang,TEOS_bonds%cang)
!      mol_atoms = (/ iSilicon,iOxygenH,iOxygenH /)
      allocate( mol_vec(3,n_atom_GAS) )
mol_vec = pr_TEOS%r
! mol_vec(:,1) = 0.0_wp
! mol_vec(:,2) = (/ 0.0_wp, 0.0_wp, ASiO /)
! mol_vec(:,3) = (/ 0.0_wp, -(2.0_wp*((1.0_wp/3.0_wp)*ASiO))*sqrt(2.0_wp), (1.0_wp/3.0_wp)*ASiO /)
!      if (n_atom_GAS > 0) call tetra_coords5(ASiO, mol_vec)

print '(3f12.6)',mol_vec
print *,mol_atoms
print '(a)',atom_name(mol_atoms)
!     Initialize the gas molecules
      do i = 1,n_GAS
         GAS_charge(i,:) = pr_TEOS%q(i)
         GAS_atom(i,:) = mol_atoms(1:n_atom_gas)
         do j = 1, n_atom_gas
            GAS_mass(i,j,:) = amass(mol_atoms(j))
         end do
      end do
      GAS_massi = 1.0_wp/GAS_mass

      rad_GAS = sigLJ_2(mol_atoms)



      allocate(mol_atoms_h2o(n_atom_h2o))
      allocate(h2o_atom(1:n_h2o,n_atom_h2o))
      allocate(h2o_charge(n_h2o,n_atom_h2o))
      allocate(h2o_mass(1:n_h2o,n_atom_h2o,3),h2o_massi(1:n_h2o,n_atom_h2o,3))
      allocate(rad_h2o(n_atom_h2o))
      if (n_atom_h2o == 1) mol_atoms_h2o = (/ iOw /)
      allocate( mol_vec_h2o(3,n_atom_h2o) )
      if (n_atom_h2o == 1) mol_vec_h2o(:,1) = 0.0_wp
! mol_vec_h2o(:,2) = (/ 0.0_wp, 0.0_wp, ??? /)
! mol_vec_h2o(:,3) = (/ 0.0_wp, ???, ??? /)
!if (n_atom_h2o > 0) call h2o_coords(???,???)
print *,'H2O'
print '(3f12.6)',mol_vec_h2o
print *,mol_atoms_h2o
print '(a)',atom_name(mol_atoms_h2o)
!     Initialize the h2o molecules
      do i = 1,n_h2o
         !h2o_charge(i,:) =
         h2o_atom(i,:) = mol_atoms_h2o(1:n_atom_h2o)
         do j = 1, n_atom_h2o
            h2o_mass(i,j,:) = amass(mol_atoms_h2o(j))
         end do
      end do
if(n_atom_h2o==1) h2o_mass = mOx+2*mH
      h2o_massi = 1.0_wp/h2o_mass

      rad_h2o = sigLJ_2(mol_atoms_h2o)

!
!--------------------------------NB-----change nat = 3 if electrosttics used ----------
!     nat = 3
!     nat = 2
!
      allocate(h2o_xyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(r_h2o_uc(n_h2o,n_atom_h2o,3))
      allocate(dr_h2o(n_h2o,n_atom_h2o,3))
      allocate(rcm_h2o_uc(n_h2o,3))
      allocate(vcm_h2o(n_h2o,3))
      allocate(image_h2o_xyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(h2o_fxyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
      allocate(h2o_vxyz(1:n_h2o,1:n_atom_h2o,1:NDIM))
!
      allocate(GAS_xyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(r_GAS_uc(n_GAS,n_atom_gas,3))
      allocate(dr_GAS(n_GAS,n_atom_gas,3))
      allocate(rcm_GAS_uc(n_GAS,3))
      allocate(vcm_GAS(n_GAS,3))
      allocate(image_GAS_xyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(GAS_fxyz(1:n_GAS,1:n_atom_gas,1:NDIM))
      allocate(GAS_vxyz(1:n_GAS,1:n_atom_gas,1:NDIM))
!
      allocate(Ox_xyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(r_O2_uc(n_O2,n_atom_o2,3))
      allocate(dr_O2(n_O2,n_atom_o2,3))
      allocate(rcm_O2_uc(n_O2,3))
      allocate(vcm_O2(n_O2,3))
      allocate(image_Ox_xyz(1:n_O2,1:NDIM,1:NDIM))
      allocate(Ox_fxyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(Ox_vxyz(1:n_O2,1:n_atom_o2,1:NDIM))
      allocate(Ox_atom(1:n_O2,3))
      allocate(Ox_gas_charge(n_O2,3))
      allocate(Ox_mass(1:n_O2,n_atom_o2,3))
      allocate(Ox_massi(1:n_O2,n_atom_o2,3))
!
      allocate(N2_xyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(r_N2_uc(n_N2,n_atom_n2,3))
      allocate(dr_N2(n_N2,n_atom_n2,3))
      allocate(rcm_N2_uc(n_N2,3))
      allocate(vcm_N2(n_N2,3))
      allocate(image_N2_xyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_fxyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_vxyz(1:n_N2,1:n_atom_n2,1:NDIM))
      allocate(N2_atom(1:n_N2,3))
      allocate(N2_gas_charge(n_N2,3))
      allocate(N2_mass(1:n_N2,n_atom_n2,3))
      allocate(N2_massi(1:n_N2,n_atom_n2,3))
!
      allocate(CO2_xyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(r_CO2_uc(n_CO2,n_atom_co2,3))
      allocate(dr_CO2(n_CO2,n_atom_co2,3))
      allocate(rcm_CO2_uc(n_CO2,3))
      allocate(vcm_CO2(n_CO2,3))
      allocate(image_CO2_xyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_fxyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_vxyz(1:n_CO2,n_atom_co2,1:NDIM))
      allocate(CO2_atom(1:n_CO2,n_atom_co2))
      allocate(CO2_gas_charge(n_CO2,n_atom_co2))
      allocate(CO2_mass(1:n_CO2,n_atom_co2,3))
      allocate(CO2_massi(1:n_CO2,n_atom_co2,3))

      allocate( ULJEL(max(n_o2,n_n2,n_co2,n_gas,n_h2o)) )
      allocate( e_bond_CO2(n_co2),e_angle_CO2(n_co2),e_bond_gas(n_gas),e_angle_gas(n_gas) )



      select case(ftyp)
      case('xyz')
         read (io_config,*) natom
         read (io_config,*) boxl
         io_topol = 15
         call get_free_file_unit(io_topol)
         open(unit=io_topol,file=trim(coordfile(1:nc - 3))//'con')
      case('pdb')
         read (io_config,*) c6,natom  ;print *,'natom = ',natom
         read (io_config,*) c6,boxl   ;print '(a,6f14.8)',' boxl = ',boxl
         ! io_topol = io_config
      case default
         print *,"'",ftyp,"'"
         stop 'unknown file type'
      end select

!
      natom_max = natom
      boxl = boxl/Angstrom
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      volbox = boxl(1)*boxl(2)*boxl(3)*(Angstrom**3) ! Box volume in A^3
!rcutel = boxl2
!rcutel2 = boxl2**2
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(fxyz(1:natom_max,3),vxyz(1:natom_max,3))
      allocate(atom(0:natom_max),latom(natom_max),atom2(natom_max))
      atom = 0
      allocate(charge(natom))
      allocate(proximity(natom_max,4))
      proximity = 0
      allocate(mass(natom_max,3),massi(natom_max,3))


      select case(ftyp)
      case('xyz')
         do i = 1,natom
            read (io_config,*) ctmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
         end do
         close(io_config)
      case('pdb')
         do i = 1,natom
            read (io_config,'(a6)',advance = 'no') c6
            read (io_config,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
            atom(i) = name2atom(trim(ctmp))
!print '(i6,3f12.6,i6)',i,rxyz(i,:),atom(i)
         end do
      end select
      rxyz(1:natom,:) = rxyz(1:natom,:)/Angstrom

      do i = 1,natom
         mass(i,:) = amass(atom(i))
      end do
      massi(1:natom,:) = 1.0_wp/mass(1:natom,:)


      select case(ftyp)
      case('xyz')
         call read_conn(io_topol)
!         read(io_topol,*) itmp
!         if(itmp /= natom) stop 'itmp /= natom'
!         do i = 1,natom
!            read(io_topol,'(a32)') ctmp
!            do j = 1,ncmax(atom(i))
!               c5 = ctmp(6 + 5*j + 1:6 + 5*j + 5)
!               read( unit=c5,fmt=* ) proximity(i,j)
!            end do
!         end do
         close(io_topol)
      case('pdb')
         do i = 1,natom
            read(io_config,'(a32)') ctmp
            do j = 1,ncmax(atom(i))
               c5 = ctmp(6 + 5*j + 1 : 6 + 5*j + 5)
               read( unit=c5,fmt=* ) proximity(i,j)
            end do
         end do
         close(io_config)
      end select

      call assign_charge(1,natom)

      write(*,*) 'sum charge = ',sum(charge(1:natom))
      write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
allocate(listSph(natom))
allocate(listALL(natom))
allocate(alst(natom_max))
listALL = (/ (i,i=1,natom) /)
      n_si = count(atom(1:natom) == iSilicon)
      n_ob = count(atom(1:natom) == iOxygen)
      n_oh = count(atom(1:natom) == iOxygenH)
      write(*,*) 'ntot = ',n_si + n_ob + n_oh
      write(*,*) '# Si = ',n_si
      write(*,*) '# O  = ',n_ob
      write(*,*) '# OH = ',n_oh
      write(*,*) 'number % = ',100.0*(n_si + n_ob + n_oh)/24000.0_wp
      write(*,*) 'density % = ',100.0*(n_si*mSi + n_ob*mOx + n_oh*mOh)/(8000*mSi + 16000*mOx)
      fdens_wt = frac_dens_wt()
      write(*,*) 'density % = ',100.0*fdens_wt
      write(*,*) 'N_Si = ',n_si/volbox,' Angstrom^-3'
      write(*,*) 'N_O  = ',n_ob/volbox,' Angstrom^-3'
      write(*,*) 'N_OH = ',n_oh/volbox,' Angstrom^-3'
      dens_gcm3 = (n_si*mSi + n_ob*mOx + n_oh*mOh)*(1e24_wp/(volbox*NAvo))
      write(*,*) 'density = ',dens_gcm3,' g/cm^3'

!allocate(listOH(count(atom == iOxygenH)))
!n_OH = 0
!do i = 1,natom
!   if (atom(i) == iOxygenH) then
!      n_OH = n_OH + 1
!      listOH(n_OH) = i
!   end if
!end do
!print *,'n_OH',n_OH
!STOP
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
!


!      ! translate the box by -L_z/2
!      rxyz(1:natom,3) = rxyz(1:natom,3) - boxl2(3)   ;print *,'boxl(3) ',boxl(3)
!      boxl(3) = 2*boxl(3)         ;print *,'boxl(3) ',boxl(3)
!      boxl2(3) = boxl(3)/2.0_wp   ;print *,'boxl2(3) ',boxl2(3)
!      boxli(3) = 1.0_wp/boxl(3)   ;print *,'boxli(3) ',boxli(3)
!
!call new_frame_file(imve,'frame',10000)
!call write_box_frame(imve,n_CO2,print_all = .true.)
call check_proximity(1,natom,ok,i)
if (.NOT.ok) then
   write(*,*) 'proximity array is NOT consistent'
   write(*,*) 'check_proximity ',ok,i
   stop
end if
call check_proximity2(1,natom,atom(1:natom),proximity,ok,i)
if (.NOT.ok) then
  write(*,*) 'proximity array is NOT consistent'
  write(*,*) 'check_proximity2 ',ok,i
  call write_atom_info(6,i)
  stop
end if

!
!      dLflexible = 10.0_wp/Angstrom
!      dLfixed = 5.0_wp/Angstrom
!      dlkeep = dLflexible + dLfixed
!      i = 0
!      do
!         i = i + 1
!         if(i > natom) exit
!         if (rxyz(i,3) < -dlkeep) then
!            call delete_atom(i)  ! ; print *,'delete ',i
!            i = i - 1
!         end if
!      end do
!
!call new_frame_file(imve,'frame',20000)
!call write_box_frame(imve,n_CO2,print_all = .true.)
!
!      nlst = 0
!      do i = 1,natom
!         if (rxyz(i,3) < -dLflexible) then
!            nlst = nlst + 1
!            alst(nlst) = i
!         end if
!      end do
!      do i = 1,nlst
!         call swap_atom(i,alst(i))   ; print *,'swap ',i,alst(i)
!      end do
!      nfixed = nlst
!      print *,'nfixed = ',nfixed
!      print *,' natom = ',natom


!      do i = 1,natom
!         if (rxyz(i,3) >= -dLflexible) then
!            do j = i+1, natom
!               if (rxyz(j,3) < -dLflexible) then
!                  call swap_atoms(i,j); print *,'swap ',i,j
!                  exit
!               end if
!            end do
!         end if
!      end do

!atom(1:nfixed) = iOw
!call new_frame_file(imve,'frame',30000)
!call write_box_frame(imve,n_CO2,print_all = .true.)

call check_proximity(1,natom,ok,i)
if (.NOT.ok) then
   write(*,*) 'proximity array is NOT consistent'
   write(*,*) 'check_proximity ',ok,i
end if
!call check_proximity2(1,natom,ok,i)
!if (.NOT.ok) then
!  write(*,*) 'proximity array is NOT consistent'
!  write(*,*) 'check_proximity2 ',ok,i
!  stop
!end if
!
!
!      do i = nfixed+1,natom
!         if (ALL(proximity(i,:)==0)) then
!         !if (count(proximity(i,:) > 0) == 0) then
!            print *,'free atom: ',i,atom_name(atom(i)),count(proximity(i,:) > 0)
!         end if
!         select case(atom(i))
!         case(iOxygen)
!            if (count(proximity(i,:) > 0) < 2) then
!               print *,i,atom_name(atom(i)),count(proximity(i,:) > 0)
!            end if
!         case(iOxygenH)
!         case(iSilicon)
!            if (count(proximity(i,:) > 0) < ncmax(atom(i))) then
!               print *,i,atom_name(atom(i)),count(proximity(i,:) > 0)
!            end if
!            nc = count( atom(proximity(i,:))==iOxygen )
!            if (nc < ncmax(atom(i))) then
!               print *,i,atom_name(atom(i)),nc,noh(i),nc+noh(i),count(proximity(i,:) > 0)
!               print '(i5,2x,a2,2x,4i6,4(1x,a2))',i,atom_name(atom(i)),proximity(i,:),atom_name(atom(proximity(i,:)))
!            end if
!         case default
!            print *,i,atom(i)
!            stop 'unknown atom type'
!         end select
!      end do
!!
!call new_frame_file(imve,'frame',30000)
!call write_box_frame(imve,n_CO2,print_all = .true.)
!call write_frame(imve,1,natom)
!! remove the free oxygens
!         j = natom
!         noxdel = 0
!         do
!            if (atom(j) == iOxygen .or. atom(j) == iOxygenH) then
!            if (count(proximity(j,:) > 0) < 1) then
!               call delete_atom(j)
!               noxdel = noxdel + 1
!            end if
!            end if
!            j = j - 1
!            if (j <= 0) exit
!         end do
!! convert singly bonded oxygens to OH
!         do i = 1,natom
!            if (atom(i) == iOxygen) then
!            if (count(proximity(i,:) > 0) < 2) then
!               atom(i) = iOxygenH
!            end if
!            end if
!         end do

!stop

      call INIT_NLIST(boxl(1),boxl(2),boxl(3),del_nlist,natom_max)
      call NEW_NLIST(1,natom)
      call Init_Verlet_List(rv0=5.0_wp/Angstrom,rcut0=2.6_wp/Angstrom)
      call new_vlist


!
      do i = 1,n_O2
        Ox_gas_charge(i,1) = q_O_O2
        Ox_gas_charge(i,2) = q_O_O2
        Ox_gas_charge(i,3) = -2.0_wp*q_O_O2
        Ox_atom(i,1) = iO_O2
        Ox_atom(i,2) = iO_O2
        Ox_mass(i,:,:) = mOx
      end do
      Ox_massi = 1.0_wp/Ox_mass
      do i = 1,n_N2
        N2_gas_charge(i,1) = q_n_N2
        N2_gas_charge(i,2) = q_n_N2
        N2_gas_charge(i,3) = -2.0_wp*q_n_N2
        N2_atom(i,1) = in_N2
        N2_atom(i,2) = in_N2
        N2_mass(i,:,:) = mN
      end do
      N2_massi = 1.0_wp/N2_mass
      do i = 1,n_CO2
        CO2_gas_charge(i,1) = q_C_CO2
        CO2_gas_charge(i,2) = q_O_CO2
        CO2_gas_charge(i,3) = q_O_CO2
        CO2_atom(i,1) = iC_CO2
        CO2_atom(i,2) = iO_CO2
        CO2_atom(i,3) = iO_CO2
        CO2_mass(i,1,:) = mC
        CO2_mass(i,2:3,:) = mOx
      end do
      CO2_massi = 1.0_wp/CO2_mass

! Insertion
rl = (/ -boxl2(1),-boxl2(2),0.0_wp /)
ru =  (/ boxl2(1), boxl2(2),boxl2(3)-10.0_wp/Angstrom /)

      call insert_mol_GAS(n_h2o,n_atom_h2o,mol_vec_h2o,h2o_xyz,100000,rad_h2o,rl,ru,fover,check_self_overlap) ; print *,'n_h2o = ',n
      call insert_mol_GAS(n_h2o,n_atom_h2o,mol_vec_h2o,h2o_xyz,100000,rad_h2o,rl,ru,fover,check_self_overlap) ; print *,'n_h2o = ',n

      call insert_mol_GAS(n_GAS,n_atom_gas,mol_vec,GAS_xyz,100000,rad_GAS,rl,ru,fover,check_self_overlap) ; print *,'n_gas = ',n_gas
      if (allocated(mol_vec)) deallocate(mol_vec)
      if (.not.allocated(mol_vec)) allocate(mol_vec(3,3))

      if (n_O2 > 0) then
      mol_vec(1,1) =  AOO*0.5_wp
      mol_vec(1,2) = -AOO*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_O2,n_atom_o2,mol_vec,Ox_xyz,1000000,rad_Ox,rl,ru,fover,check_self_overlap); print *,'n_o2 = ',n_o2
      end if
      if (n_N2 > 0) then
      mol_vec(1,1) =  ANN*0.5_wp
      mol_vec(1,2) = -ANN*0.5_wp
      mol_vec(1,3) = 0.0_wp
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_N2,n_atom_n2,mol_vec,N2_xyz,1000000,rad_N2,rl,ru,fover,check_self_overlap); print *,'n_n2 = ',n_n2
      end if
      if (n_CO2 > 0) then
      mol_vec(1,1) = 0.0_wp
      mol_vec(1,2) =  ACO2
      mol_vec(1,3) = -ACO2
      mol_vec(2,:) = 0.0_wp
      mol_vec(3,:) = 0.0_wp
      call insert_mol_GAS(n_CO2,n_atom_co2,mol_vec,CO2_xyz,1000000,rad_CO2,rl,ru,fover,check_self_overlap); print *,'n_co2 = ',n_co2
      end if
!
      if (restart) then
         call read_box_frame(io_restart_xyz,nmol = n_GAS,print_all = .true.)
         call read_vel(io_restart_v,nmol = n_GAS,print_all = .true.)
      end if
!
! Set up the uncorrected coords
      do i = 1, n_h2o
         r_h2o_uc(i,1,:) = h2o_xyz(i,1,:)
         do j = 2,n_atom_h2o
            rr = h2o_xyz(i,j,:) - h2o_xyz(i,1,:)
            call pbc(rr)
            r_h2o_uc(i,j,:) =  r_h2o_uc(i,1,:) + rr
         end do
      end do

      do i = 1, n_GAS
         r_GAS_uc(i,1,:) = GAS_xyz(i,1,:)
         do j = 2,n_atom_GAS
            rr = GAS_xyz(i,j,:) - GAS_xyz(i,1,:)
            call pbc(rr)
            r_GAS_uc(i,j,:) =  r_GAS_uc(i,1,:) + rr
         end do
      end do
      do i = 1, n_O2
         r_O2_uc(i,1,:) = Ox_xyz(i,1,:)
         rr = Ox_xyz(i,2,:) - Ox_xyz(i,1,:)
         call pbc(rr)
         r_O2_uc(i,2,:) = r_O2_uc(i,1,:) + rr
      end do
      do i = 1, n_N2
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
         call write_box_frame(imve,n_CO2,print_all = .true.)
      end if
!========================================================================
!
      call LJ_INIT

!      call comvel(natom,T_kelvin,vxyz,mass)
!      call comvel(n_h2o,n_atom_h2o,T_kelvin,h2o_vxyz,h2o_mass)
!      call comvel(n_GAS,n_atom_GAS,T_kelvin,GAS_vxyz,GAS_mass)
!!GAS_vxyz(:,:,3) = -0.01
!      call comvel(n_O2,2,T_kelvin,Ox_vxyz,Ox_mass)
!      call comvel(n_N2,2,T_kelvin,N2_vxyz,N2_mass)
!      call comvel(n_CO2,3,T_kelvin,CO2_vxyz,CO2_mass)
!!
!      print '(a/,3g18.8)','net GAS momentum', &
!         sum(GAS_vxyz(:,:,1)*GAS_mass(:,:,1)),&
!         sum(GAS_vxyz(:,:,2)*GAS_mass(:,:,2)),&
!         sum(GAS_vxyz(:,:,3)*GAS_mass(:,:,3))
!      SUMX = 0.0_wp
!      SUMY = 0.0_wp
!      SUMZ = 0.0_wp
!      do j = 1,n_atom_GAS
!      do i = 1,n_GAS
!         SUMX = SUMX + GAS_vxyz(i,j,1)*GAS_mass(i,j,1)
!         SUMY = SUMY + GAS_vxyz(i,j,2)*GAS_mass(i,j,2)
!         SUMZ = SUMZ + GAS_vxyz(i,j,3)*GAS_mass(i,j,3)
!      end do
!      end do
!      print '(3g18.8)',SUMX,SUMY,SUMZ
!      print *
!
!      print '(a/,3g18.8)','net H2O momentum', &
!         sum(h2o_vxyz(:,:,1)*h2o_mass(:,:,1)),&
!         sum(h2o_vxyz(:,:,2)*h2o_mass(:,:,2)),&
!         sum(h2o_vxyz(:,:,3)*h2o_mass(:,:,3))
!
!      print '(a/,3g18.8)','net CO2 momentum', &
!         sum(CO2_vxyz(:,:,1)*CO2_mass(:,:,1)),&
!         sum(CO2_vxyz(:,:,2)*CO2_mass(:,:,2)),&
!         sum(CO2_vxyz(:,:,3)*CO2_mass(:,:,3))
!!print *,'CO2_mass(:,1,1:NDIM))',sum(CO2_mass(:,1,1:NDIM))/(n_CO2*3)
!!print *,'CO2_mass(:,2,1:NDIM))',sum(CO2_mass(:,2,1:NDIM))/(n_CO2*3)
!!print *,'CO2_mass(:,3,1:NDIM))',sum(CO2_mass(:,3,1:NDIM))/(n_CO2*3)
!
!      SUMX = 0.0_wp
!      SUMY = 0.0_wp
!      SUMZ = 0.0_wp
!      do j = 1,3
!      do i = 1,n_CO2
!         SUMX = SUMX + CO2_vxyz(i,j,1)*CO2_mass(i,j,1)
!         SUMY = SUMY + CO2_vxyz(i,j,2)*CO2_mass(i,j,2)
!         SUMZ = SUMZ + CO2_vxyz(i,j,3)*CO2_mass(i,j,3)
!      end do
!      end do
!      print '(3g18.8)',SUMX,SUMY,SUMZ
!      print *
!
!      print '(a/,3g18.8)','net solid momentum', &
!         sum(vxyz(:,1)*mass(:,1)),&
!         sum(vxyz(:,2)*mass(:,2)),&
!         sum(vxyz(:,3)*mass(:,3))
!      print '(a/,3g18.8)','net solid momentum', &
!         dot_product(vxyz(:,1),mass(:,1)),&
!         dot_product(vxyz(:,2),mass(:,2)),&
!         dot_product(vxyz(:,3),mass(:,3))
!
!
!!      call voidage_calc(probe_O2,-boxl2,boxl2,100000,voidage(1))
!!      call voidage_calc(probe_N2,-boxl2,boxl2,100000,voidage(2))
!!      call voidage_calc(probe_CO2,-boxl2,boxl2,100000,voidage(3))
!     ! call voidage_calc(pr_teos,-boxl2,boxl2,100000,voidage(3))
!call start_clock(t(1))
!      !call voidage_calc(pr_teos,-boxl2,boxl2,10000,voidage(4))
!      call voidage_calc(probe_SiO4,-boxl2,boxl2,10000,voidage(4))
!      !call voidage_calc(probe_tip3p,-boxl2,boxl2,10000,voidage(4))
!call stop_clock(t(1))
!call write_statistic(t(1)%dt,6)
!print *,'         voidage = ',voidage(4)
!      !call voidage_calc(probe_tip3p,-boxl2,boxl2,100000,voidage(5))j
!  !call voidage_profile(probe_O2,-boxl2,boxl2,10000,5.0_wp/Angstrom,nbinv,voidage)
nb=int(boxl/0.5_wp)
print *,'nb = ',nb
!nx = nb(1)
!ny = nb(2)
!nz = nb(3)
!allocate(voidxyz(nb(1),nb(2),nb(3)))
allocate(HenryK(nb(1),nb(2),nb(3)))
!call start_clock(t(1))
!!call voidage_profile(pr_teos,-boxl2,boxl2,20,0.5_wp,nb,voidxyz)
!call voidage_profile(probe_SiO4,-boxl2,boxl2,5,0.5_wp,nb,voidxyz)
!!call voidage_profile(probe_tip3p,-boxl2,boxl2,125,0.5_wp,nb,voidxyz)
!call stop_clock(t(1))
!call write_statistic(t(1)%dt,6)
!print *,'nb = ',nb
!print *,'sum voidxyz/nbin = ',sum(voidxyz)/(nb(1)*nb(2)*nb(3))
!print *,' maxval(voidxyz) = ',maxval(voidxyz)
!print *,' minval(voidxyz) = ',minval(voidxyz)
!open(unit=93, file='voidxyz.out',status='unknown')
!write(93,*)nb(1:3)
!do iz = 1,nb(3)
!do iy = 1,nb(2)
!do ix = 1,nb(1)
!   write(93,'((1x,f0.8))',advance='no') voidxyz(ix,iy,iz)
!end do
!write(93,*)
!end do
!end do
!close(93)


rcutH = 10.0_wp/Angstrom
ntrial = 100000
HenryK = 0.0_wp
dl = 5.0_wp/Angstrom
rl = -boxl2
ru = (/ -boxl2(1)+dl,-boxl2(2)+dl,-boxl2(3)+dl /)
!print *,'# probe_sio4q 1'
!call Henrys_law_calc2(probe_sio4q,rl,ru,ntrial,0.8_wp,rcutH,Uljel_1,Khenry,facc)
!print '(4g16.8)',Khenry,Uljel_1,facc
!print *
!print *
!
!print *,'# probe_sio4q 2'
!call Henrys_law_calc2(probe_sio4q,rl,ru,ntrial,0.8_wp,rcutH,Uljel_1,Khenry,facc)
!print '(4g16.8)',Khenry,Uljel_1,facc
!print *
print *
print *,'# probe_sioh4 1'
call Henrys_law_calc2(probe_sioh4,rl,ru,ntrial,0.8_wp,rcutH,Uljel_1,Khenry,facc)
print '(4g16.8)',Khenry,Uljel_1,facc
print *
print *

rl = (/ -boxl2(1)+dl,-boxl2(2),-boxl2(3) /)
ru = (/ -boxl2(1)+2*dl,-boxl2(2)+dl,-boxl2(3)+dl /)

print *,'# probe_sioh4 2'
call Henrys_law_calc2(probe_sioh4,rl,ru,ntrial,0.8_wp,rcutH,Uljel_1,Khenry,facc)
print '(4g16.8)',Khenry,Uljel_1,facc
print *
print *



stop



rcutH = 10.0_wp/Angstrom
ntrial = 100
open(unit=94, file='henryk.out',status='unknown')

do k = 1,12
print *,'rcutH = ',rcutH

HenryK = 0.0_wp
dl = 5.0_wp/Angstrom
call start_clock(t(1))
rl = -boxl2
ru = (/ -boxl2(1)+2*dl,-boxl2(2)+dl,-boxl2(3)+dl /)
   call Henrys_law_profile(probe_sioh4,rl,ru,ntrial,0.8_wp,rcutH,dl,nb,HenryK)
call stop_clock(t(1))
call write_statistic(t(1)%dt,6)
print *,'nb = ',nb
print *,'sum HenryK/nbin = ',sum(HenryK)/(nb(1)*nb(2)*nb(3))
print *,' maxval(HenryK) = ',maxval(HenryK)
print *,' minval(HenryK) = ',minval(HenryK)
!write(94,*)nb(1:3)
write(94,'(3(1x,i9))')nb
do iz = 1,nb(3)
do iy = 1,nb(2)
do ix = 1,nb(1)
   write(94,'((1x,f0.8))',advance='no') HenryK(ix,iy,iz)
end do
write(94,*)
end do
end do

rcutH = rcutH+dl
ntrial = ntrial * 5
end do

stop




!open(unit=93, file='voidxyz.out',status='unknown')
!write(93,*)nb(1:3)
!do iz = 1,nb(3)
!do iy = 1,nb(2)
!do ix = 1,nb(1)
!   write(93,'(3(i4,1x),f0.8)') ix,iy,iz,voidxyz(ix,iy,iz)
!end do
!end do
!end do
!close(93)

call start_clock(t(1))

call hk_initialize((nb(1)*nb(2)*nb(3))/2)
allocate(Lattice(nb(1),nb(2),nb(3)))
allocate(Lb(nb(1),nb(2),nb(3)),lab((nb(1)*nb(2)*nb(3))/2))


where(voidxyz >= 0.035)
   Lattice=1
elsewhere
   Lattice=0
end where

!open(unit=93, file='voidxyz.out',status='unknown')
!write(93,*)nb(1:3)
!do iz = 1,nb(3)
!do iy = 1,nb(2)
!do ix = 1,nb(1)
!   write(93,'((1x,i1))',advance='no') Lattice(ix,iy,iz)
!end do
!write(93,*)
!end do
!!write(93,*)
!end do
!close(93)

call print_mat(6,Lattice); write(6,'(50("-"))')

call hk_lattice_3d(Lattice,periodic=(/.true.,.true.,.false./),nclusters=nclusters)
!call printing("Labeled_clusters.txt")
!call print_mat(6); write(6,'(50("-"))')


do iz = 1,nb(3)
write(*,*) 'iz = ',iz
do iy = 1,nb(2)
do ix = 1, nb(1)
   if (Lattice(ix,iy,iz) /= 0) then
      write(*,"((1x,i3))",advance='no') Lattice(ix,iy,iz)
   else
      write(*,"(a4)",advance='no') '  _ '
   end if
end do
write(*,*)
end do
write(*,*)
end do

! percolation check.
print *,'percolation cluster x = ',percolation_check(Lattice,dim_xyz=1)
print *,'percolation cluster y = ',percolation_check(Lattice,dim_xyz=2)
print *,'percolation cluster z = ',percolation_check(Lattice,dim_xyz=3)

call stop_clock(t(1))
call write_statistic(t(1)%dt,6)

Lb = .false.
nlab = 0
iz = nb(3)
do iy = 1,nb(2)
do ix = 1,nb(1)
   k = Lattice(ix,iy,iz)
   if (k==0)cycle
   if(.not.in_list(k,lab,nlab)) then
      nlab = nlab + 1
      lab(nlab) = k
   end if
end do
end do
print *,'nlab = ',nlab
print *,'lab = ',lab(1:nlab)
print *
do iz = nb(3),1,-1
do iy = 1,nb(2)
do ix = 1,nb(1)
   k = Lattice(ix,iy,iz)
   if (in_list(k,lab,nlab))then
      lb(ix,iy,iz)=.true.
   end if
end do
end do
print *,nb(3)-iz+1,count(lb(:,:,iz))
end do






!print *,'nbin = ',nbinv
print '(i6,8g18.8)',natom,fdens_wt,dens_gcm3,voidage(3:5)
open(unit=193, file='voidage.out',status='unknown',position='append')
write(193,'(i6,8g18.8)')natom,fdens_wt,dens_gcm3,voidage(3:5)

stop





! Rescale each gas to the set temperature
      Ek = kinetic_energy(natom,vxyz,mass)
      Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
      Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
      Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
      Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
      Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

      if (n_h2o > 0) then
         Ndf = 3*(n_atom_h2o*n_h2o)
         T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [h2o] = ',T_kinetic
         h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [h2o] = ',T_kinetic, 'rescaled'
      end if

      if (n_GAS > 0) then
         Ndf = 3*(n_atom_gas*n_GAS)
         T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Gas] = ',T_kinetic
         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [Gas] = ',T_kinetic, 'rescaled'
      end if

      if (n_co2 > 0) then
         Ndf = 3*(3*n_CO2)
         T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [CO2] = ',T_kinetic
         CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)
         T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [CO2] = ',T_kinetic, 'rescaled'
      end if

      if (n_o2 > 0) then
         Ndf = 3*(2*n_O2)
         T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [O2] = ',T_kinetic
         Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [O2] = ',T_kinetic, 'rescaled'
      end if

      if (n_n2 > 0) then
         Ndf = 3*(2*n_N2)
         T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [N2] = ',T_kinetic
         N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
         print *,'T_kinetic [N2] = ',T_kinetic, 'rescaled'
      end if

      if (natom > 0) then
      Ndf = 3*(natom)
      T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
      print *,'T_kinetic [Solid] = ',T_kinetic
      vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
      Ek = kinetic_energy(natom,vxyz,mass)
      T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
      print *,'T_kinetic [Solid] = ',T_kinetic, 'rescaled'
      end if

      if (make_movie) then ! write out first frame of velocities
         iframe = 0
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_CO2,print_all = .false.)
      end if

      call set_bond_angle_lists()
      call new_tcf(n_h2o,ntau_msd,nsamp,msd_h2o,'msd_h2o.out')
      call new_tcf(n_h2o,ntau_vcf,nsamp,vcf_h2o,'vcf_h2o.out')
      call new_tcf(n_GAS,ntau_msd,nsamp,msd_gas,'msd_GAS.out')
      call new_tcf(n_GAS,ntau_vcf,nsamp,vcf_gas,'vcf_GAS.out')
      call new_tcf(n_O2,ntau_msd,nsamp,msd_o2,'msd_O2.out')
      call new_tcf(n_O2,ntau_vcf,nsamp,vcf_o2,'vcf_O2.out')
      call new_tcf(n_N2,ntau_msd,nsamp,msd_n2,'msd_N2.out')
      call new_tcf(n_N2,ntau_vcf,nsamp,vcf_n2,'vcf_N2.out')
      call new_tcf(n_CO2,ntau_msd,nsamp,msd_co2,'msd_CO2.out')
      call new_tcf(n_CO2,ntau_vcf,nsamp,vcf_co2,'vcf_CO2.out')

      time = 0.0_wp

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      h2o_fxyz = 0.0_wp
      GAS_fxyz = 0.0_wp
      energy = 0.0_wp
!
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL) ; energy = energy + ULJEL
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_h2o,n_atom_h2o,h2o_atom,h2o_xyz,h2o_fxyz,ULJEL); print*,'ULJ h2o  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_h2o))
do i = 1,n_h2o
 print *,i,ULJEL(i)
end do
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL); print*,'ULJ GAS  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_GAS))
do i = 1,n_gas
 print *,i,ULJEL(i)
end do

      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL); print*,'ULJ O2  ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_O2))
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL); print*,'ULJ N2 ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_N2))
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL); print*,'CO2 GAS ',sum(ULJEL)
      energy = energy + sum(ULJEL(1:n_CO2))
      call force_keating(ukeating)
      call repul_force2(1,natom,repulenergy)


      call O2_stretching(e_bond_o2) ; energy = energy + e_bond_o2
      call N2_stretching(e_bond_n2) ; energy = energy + e_bond_n2
      call CO2_stretching(e_bond_co2) ; energy = energy + sum(e_bond_co2(1:n_co2))

      do i = 1, n_GAS
         call force_keating_bond2(TEOS_bonds%nbond,TEOS_bonds%ibond,TEOS_bonds%kbond,TEOS_bonds%abond,&
                                  GAS_xyz(i,:,:),GAS_fxyz(i,:,:),e_bond_teos)
         energy = energy + e_bond_teos
         call force_keating_angle2(TEOS_bonds%nang,TEOS_bonds%iang,TEOS_bonds%kang,TEOS_bonds%cang,&
                                   GAS_xyz(i,:,:),GAS_fxyz(i,:,:),e_angle_teos)
         energy = energy + e_angle_teos
      end do


      call CO2_angle_bending(e_angle_co2)
      energy = energy + sum(e_angle_co2(1:n_co2))
!      if( n_atom_h2o > 1 )then
!         call GAS_stretching(n_h2o,n_atom_h2o,h2o_xyz,A_h2o,k_ohw,e_bond_h2o)
!         energy = energy + sum(e_bond_h2o(1:n_h2o))
!         call GAS_angle_bending(n_h2o,n_atom_h2o,h2o_xyz,ctheta = ctheta_h2o,ktheta = K_h2o,en=e_angle_h2o)
!         energy = energy + sum(e_angle_h2o(1:n_h2o))
!      end if

print *,'e_bond_o2 ',e_bond_o2
print *,'e_bond_n2 ',e_bond_n2
print *,'e_bond_teos ',e_bond_teos
print *,'e_angle_teos ',e_angle_teos
!do i = 1,n_gas
! print *,i,e_bond_gas(i),e_angle_gas(i)
!end do



if(natom > 0)then
call list_atoms_in_sphere(rxyz(natom,:),15.0_wp/Angstrom,nSph,listSph)
print *,'No. of atoms in sphere is ',nSph
end if

time_loop: do i = 1,10
print '(//)'
call start_clock(t(1))
energy = energy4(1,natom)
call stop_clock(t(1))
print *,'            energy4 = ',energy

call start_clock(t(2))
energy = keating_energy_List(natom,listALL)
call stop_clock(t(2))
print *,'keating_energy_List = ',energy


call start_clock(t(3))
call energy_keating(energy)
call stop_clock(t(3))
print *,'     energy_keating = ',energy


call start_clock(t(4))
call force_keating(ukeating)
call stop_clock(t(4))
print *,'     energy_keating = ',ukeating


print '(/)'
call start_clock(t(5))
repulenergy = repul_energy(1,natom)
call stop_clock(t(5))
print *,'     repul_energy = ',repulenergy


call start_clock(t(6))
repulenergy = repul_energy_VL(1,natom)
call stop_clock(t(6))
print *,'  repul_energy_VL = ',repulenergy

call start_clock(t(7))
repulenergy = repul_energy_NL(1,natom)
call stop_clock(t(7))
print *,'  repul_energy_NL = ',repulenergy


call start_clock(t(8))
call repul_force(1,natom,repulenergy)
call stop_clock(t(8))
print *,'     repul_energy = ',repulenergy

call start_clock(t(9))
call repul_force2(1,natom,repulenergy)
call stop_clock(t(9))
print *,'     repul_energy = ',repulenergy

if (natom > 0) then

print '(//)'
fxyz = 0.0_wp
call start_clock(t(10))
call keating_force_atom(1,ukeating,fxyz(1,:))
call stop_clock(t(10))
print '(3(a,g18.8))','keating_force_atom(1),ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(11))
call force4(1,1,ukeating)
call stop_clock(t(11))
print '(3(a,g18.8))','          force4(1,1),ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(12))
call force_keating(ukeating)
call stop_clock(t(12))
print '(3(a,g18.8))','        force_keating,ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(13))
call keating_force_List(natom,listALL,ukeating)
call stop_clock(t(13))
print '(3(a,g18.8))','   keating_force_List,ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(14))
call force4(1,natom,ukeating)
call stop_clock(t(14))
print '(3(a,g18.8))','      force4(1,natom),ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(15))
call energy_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call energy_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(15))
print '(3(a,g18.8))','     energy_keating_bond/angle = ',Ubond+Uang

fxyz = 0.0_wp
call start_clock(t(16))
call force_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call force_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(16))
print '(3(a,g18.8))','      force_keating_bond/angle = ',Ubond+Uang,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(17))
energy = keating_energy_List(nSph,listSph)
call stop_clock(t(17))
print '(a,g18.8)','  (Sphere) keating_energy_List = ',energy

call start_clock(t(18))
call list_bond_angle(nSph,listSph,nbondtot,ibond,kbond,abond,nang,iang,kang,ctheta)
call stop_clock(t(18))

call start_clock(t(19))
call list_bond_angle_sort(nSph,listSph,nbondtot,ibond,kbond,abond,nang,iang,kang,ctheta)
call stop_clock(t(19))

fxyz = 0.0_wp
call start_clock(t(20))
call force_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call force_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(20))
print '(3(a,g18.8))','      force_keating_bond/angle = ',Ubond+Uang,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(21))
call energy_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call energy_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(21))
print '(3(a,g18.8))','     energy_keating_bond/angle = ',Ubond+Uang,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(22))
call keating_force_List(nSph,listSph,ukeating)
call stop_clock(t(22))
print '(3(a,g18.8))','   keating_force_List,ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(23))
call atom_bond_angle(1,nbondtot,ibond,kbond,abond,nang,iang,kang,ctheta)
call stop_clock(t(23))

fxyz = 0.0_wp
call start_clock(t(24))
call keating_force_List(1,(/1/),ukeating)
call stop_clock(t(24))
print '(3(a,g18.8))','   keating_force_List,ukeating = ',ukeating,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

fxyz = 0.0_wp
call start_clock(t(25))
ukeating = keating_energy_List(1,(/1/))
call stop_clock(t(25))
print '(3(a,g18.8))','  keating_energy_List,ukeating = ',ukeating

fxyz = 0.0_wp
call start_clock(t(26))
call force_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call force_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(26))
print '(3(a,g18.8))','      force_keating_bond/angle = ',Ubond+Uang,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(27))
call energy_keating_bond(nbondtot,ibond,kbond,abond,Ubond)
call energy_keating_angle(nang,iang,kang,ctheta,Uang)
call stop_clock(t(27))
print '(3(a,g18.8))','     energy_keating_bond/angle = ',Ubond+Uang,' fxyz(1,1) = ',fxyz(1,1),' fxyz(2,1) = ',fxyz(1,1)

call start_clock(t(28))
call energy_keating(energy)
call stop_clock(t(28))
print *,'     energy_keating = ',energy

call start_clock(t(29))
energy = keating_energy_atom(1)
call stop_clock(t(29))
print *,'energy_keating_atom = ',energy

call start_clock(t(30))
call set_bond_angle_lists()
call stop_clock(t(30))

end if

end do time_loop


print '(/)'
do i = 1,30
   write(6,'(i4,2x)',advance='no') i
   call write_statistic(t(i)%dt,6)
end do

if (natom > 0) then
print '(//)'
energy = energy4(1,1)
print *,'energy4 = ',energy

energy = keating_energy_atom(1)
print *,'keating_energy_atom = ',energy

rr = (/ 0.01, 0.01, 0.01 /)
uold = energy4(1,1)
rxyz(1,:) = rxyz(1,:) + rr
unew = energy4(1,1)
print *,'del energy4 = ',unew-uold
unew = keating_energy_atom(1)
rxyz(1,:) = rxyz(1,:) - rr
uold = keating_energy_atom(1)
print *,'del keating_energy_atom = ',unew-uold
end if

!call RELAX_MC_LIST(nr=1000,nl=n_OH,list=listOH,kbT=kBoltzT)
!call RELAX_MC(nr=1000,ifirst=1,ilast=natom,kbT=kBoltzT)
!stop
!call Steepest_Descent(nstep=1000,dRmax=0.1_wp/Angstrom,ifirst=1,ilast=natom,R=Rxyz,F=Fxyz)
!call Steepest_Descent_List(nstep=1000,dRmax=0.1_wp/Angstrom,nl=n_OH,list=listOH,R=Rxyz,F=Fxyz)
!stop

!
!======================================================== main loop
      main_loop: do istep = 1, number_moves

         call cpu_time(tp(1))

         coprxyz = rxyz
         call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
         call pbc(rxyz)
!
         call gas_vverlet_a(dt,h2o_massi,h2o_xyz,h2o_vxyz,h2o_fxyz,dr_h2o)
         call pbc(h2o_xyz)
         r_h2o_uc = r_h2o_uc + dr_h2o

         call gas_vverlet_a(dt,GAS_massi,GAS_xyz,GAS_vxyz,GAS_fxyz,dr_GAS)
         call pbc(GAS_xyz)
         r_GAS_uc = r_GAS_uc + dr_GAS

         call gas_vverlet_a(dt,Ox_massi,Ox_xyz,Ox_vxyz,Ox_fxyz,dr_O2)
         call pbc(Ox_xyz)
         r_O2_uc = r_O2_uc + dr_O2

         call gas_vverlet_a(dt,N2_massi,N2_xyz,N2_vxyz,N2_fxyz,dr_N2)
         call pbc(N2_xyz)
         r_N2_uc = r_N2_uc + dr_N2

         call gas_vverlet_a(dt,CO2_massi,CO2_xyz,CO2_vxyz,CO2_fxyz,dr_CO2)
         call pbc(CO2_xyz)
         r_CO2_uc = r_CO2_uc + dr_CO2

         ! check if the neighbour list needs to be updated
         do i = 1,natom
            if (CELL(rxyz(i,1:NDIM)) /= CELL(coprxyz(i,1:NDIM))) then
               call NEW_NLIST(1,natom)
               EXIT
            end if
         end do

         fxyz = 0.0_wp
         h2o_fxyz = 0.0_wp
         GAS_fxyz = 0.0_wp
         Ox_fxyz = 0.0_wp
         N2_fxyz = 0.0_wp
         CO2_fxyz = 0.0_wp

         energy = 0.0_wp                        ; call cpu_time(tp(2)); sumtp(1) = sumtp(1) + tp(2)-tp(1)
         call force_keating(ukeating)
         call repul_force2(nfixed+1,natom,repulenergy)


         do i = 1, n_GAS
            call force_keating_bond2(TEOS_bonds%nbond,TEOS_bonds%ibond,TEOS_bonds%kbond,TEOS_bonds%abond,&
                                     GAS_xyz(i,:,:),GAS_fxyz(i,:,:),e_bond_teos)
            energy = energy + e_bond_teos
            call force_keating_angle2(TEOS_bonds%nang,TEOS_bonds%iang,TEOS_bonds%kang,TEOS_bonds%cang,&
                                      GAS_xyz(i,:,:),GAS_fxyz(i,:,:),e_angle_teos)
            energy = energy + e_angle_teos
         end do
write(999,*)istep,energy

         call O2_stretching(e_bond_o2) ; energy = energy + e_bond_o2
         call N2_stretching(e_bond_n2) ; energy = energy + e_bond_n2
         call CO2_stretching(e_bond_co2) ; energy = energy + sum(e_bond_co2(1:n_co2))


         call CO2_angle_bending(e_angle_co2)
         energy = energy + sum(e_angle_co2(1:n_co2))
!        if( n_atom_h2o > 1 )then
!           call GAS_stretching(n_h2o,n_atom_h2o,h2o_xyz,A_h2o,k_ohw,e_bond_h2o)
!           energy = energy + sum(e_bond_h2o(1:n_h2o))
!           call GAS_angle_bending(n_h2o,n_atom_h2o,h2o_xyz,ctheta = ctheta_h2o,ktheta = K_h2o,en=e_angle_h2o)
!           energy = energy + sum(e_angle_h2o(1:n_h2o))
!        end if

!print *,'istep = ',istep
!print *,'e_bond_o2 = ',e_bond_o2
!print *,'e_bond_n2 = ',e_bond_n2
!print *,'e_bond_co2 = ',sum(e_bond_co2)
!print *,'e_angle_co2 = ',sum(e_angle_co2)
!print *,'e_bond_gas = ',sum(e_bond_gas)
!print *,'e_angle_gas = ',sum(e_angle_gas)


!
!     call FORCE_ENERGY_LJ_EL_GAS(1,natom,1,n_GAS,ULJEL) ; energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL); energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL); energy = energy + ULJEL
!     call FORCE_ENERGY_LJ_EL_CO2(1,natom,1,n_CO2,ULJEL); energy = energy + ULJEL
! only LJ
!      call FORCE_ENERGY_LJ(1,natom,1,n_GAS,3,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL)
!      call FORCE_ENERGY_LJ(1,natom,1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL)
! only LJ using neighbour list
      call FORCE_ENERGY_LJ_NLIST(1,n_h2o,n_atom_h2o,h2o_atom,h2o_xyz,h2o_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_h2o))
      call FORCE_ENERGY_LJ_NLIST(1,n_GAS,n_atom_GAS,GAS_atom,GAS_xyz,GAS_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_GAS))
      call FORCE_ENERGY_LJ_NLIST(1,n_O2,2,Ox_atom,Ox_xyz,Ox_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_O2))
      call FORCE_ENERGY_LJ_NLIST(1,n_N2,2,N2_atom,N2_xyz,N2_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_N2))
      call FORCE_ENERGY_LJ_NLIST(1,n_CO2,3,CO2_atom,CO2_xyz,CO2_fxyz,ULJEL) ; energy = energy + sum(ULJEL(1:n_CO2))
! Electrostatic
      call FORCE_ENERGY_EL(1,natom,1,n_GAS,3,GAS_charge,GAS_xyz,GAS_fxyz,UEL) ;  energy = energy + UEL
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_O2,Ox_gas_charge,Ox_xyz,Ox_fxyz,ULJEL)  ; energy = energy + ULJEL
!      call FORCE_ENERGY_EL_X2_dummy(1,natom,1,n_N2,n2_gas_charge,n2_xyz,n2_fxyz,ULJEL)  ; energy = energy + ULJEL
!      call FORCE_ENERGY_EL(1,natom,1,n_CO2,3,co2_gas_charge,co2_xyz,co2_fxyz,ULJEL)  ;  energy = energy + ULJEL
!




!call FORCE_ENERGY_LJ_ATTRACTIVE(listOH,n_OH,atom,rxyz,fxyz,ULJAtt)
!call FORCE_ENERGY_LJ_LIST(n_OH,listOH,atom,rxyz,fxyz,ULJAtt)
!print *,'ULJ OH = ',ULJAtt
!energy = energy + ULJAtt





      call cpu_time(tp(3)); sumtp(2) = sumtp(2) + tp(3) - tp(2)

!
      call vverlet_b(dt,massi,vxyz,fxyz)
      call gas_vverlet_b(dt,h2o_massi,h2o_vxyz,h2o_fxyz)
      call gas_vverlet_b(dt,GAS_massi,GAS_vxyz,GAS_fxyz)
      call gas_vverlet_b(dt,Ox_massi,Ox_vxyz,Ox_fxyz)
      call gas_vverlet_b(dt,N2_massi,N2_vxyz,N2_fxyz)
      call gas_vverlet_b(dt,CO2_massi,CO2_vxyz,CO2_fxyz)

      time = time + dtreal

      if (istep <= nequib) then
         ! print *,sum(fxyz(1:natom,:),dim=1)
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

         Ek_tot = Ek + Ek_O2 + Ek_N2 + Ek_CO2 + Ek_gas + Ek_h2o

         Ndf = 3*(natom + n_atom_gas*n_GAS  + n_atom_h2o*n_h2o + 2*n_O2 + 2*n_N2 + 3*n_CO2) - 3
         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)

         if (mod(istep,nsamp) == 0) then
            write(io_energy,'(8g18.8)') time, EK_tot, energy + ukeating + repulenergy*0.5_wp, &
               EK_tot + energy + ukeating + repulenergy*0.5_wp
            write(io_temp,'(8g18.8)') time, T_kinetic
            call flush(17); call flush(18)
         end if
         !if (T_kinetic == 0.0_wp) T_kinetic = 1.0_wp

!         vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
!         h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
!         GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
!         Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
!         N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
!         CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)

!         call rem_com_momentum(natom,vxyz,mass)
!         call rem_com_momentum(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
!         call rem_com_momentum(n_GAS,n_atom_GAS,GAS_vxyz,GAS_mass)
!         call rem_com_momentum(n_O2,2,Ox_vxyz,Ox_mass)
!         call rem_com_momentum(n_N2,2,N2_vxyz,N2_mass)
!         call rem_com_momentum(n_CO2,3,CO2_vxyz,CO2_mass)

         if (n_h2o > 0) then
            Ndf = 3*(n_atom_h2o*n_h2o)
            T_kinetic = Ek_h2o*2.0_wp/(K_ev*Ndf)
            h2o_vxyz = h2o_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_GAS > 0) then
            Ndf = 3*(n_atom_gas*n_GAS)
            T_kinetic = Ek_gas*2.0_wp/(K_ev*Ndf)
            GAS_vxyz = GAS_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_co2 > 0) then
            Ndf = 3*(3*n_CO2)
            T_kinetic = Ek_CO2*2.0_wp/(K_ev*Ndf)
            CO2_vxyz = CO2_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_o2 > 0) then
            Ndf = 3*(2*n_O2)
            T_kinetic = Ek_O2*2.0_wp/(K_ev*Ndf)
            Ox_vxyz = Ox_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (n_n2 > 0) then
            Ndf = 3*(2*n_N2)
            T_kinetic = Ek_N2*2.0_wp/(K_ev*Ndf)
            N2_vxyz = N2_vxyz*sqrt(T_kelvin/T_kinetic)
         end if

         if (natom > 0) then
            Ndf = 3*(natom)
            T_kinetic = Ek*2.0_wp/(K_ev*Ndf)
            vxyz = vxyz*sqrt(T_kelvin/T_kinetic)
         end if

      end if

      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then

         if (n_h2o > 0) then
         do i = 1,n_h2o
            mtmp = 0.0_wp
            rcm_h2o_uc(i,:) = 0.0_wp
            vcm_h2o(i,:) = 0.0_wp
            do j = 1,n_atom_h2o
               rcm_h2o_uc(i,:) = rcm_h2o_uc(i,:) + r_h2o_uc(i,j,:)*amass(h2o_atom(i,j))
               vcm_h2o(i,:) = rcm_h2o_uc(i,:) + h2o_vxyz(i,j,:)*amass(h2o_atom(i,j))
               mtmp = mtmp + amass(h2o_atom(i,j))
            end do
            rcm_h2o_uc(i,:) = rcm_h2o_uc(i,:)/mtmp
            vcm_h2o(i,:) = vcm_h2o(i,:)/mtmp
         end do
         end if

         if (n_gas > 0) then
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
         end if

         do i = 1,n_O2
            rcm_O2_uc(i,:) = ( r_O2_uc(i,1,:) + r_O2_uc(i,2,:) )*0.5_wp
            vcm_O2(i,:) = ( Ox_vxyz(i,1,:) + Ox_vxyz(i,2,:) )*0.5_wp
         end do
         do i = 1,n_N2
            rcm_N2_uc(i,:) = ( r_N2_uc(i,1,:) + r_N2_uc(i,2,:) )*0.5_wp
            vcm_N2(i,:) = ( N2_vxyz(i,1,:) + N2_vxyz(i,2,:) )*0.5_wp
         end do
         do i = 1,n_CO2
            rcm_CO2_uc(i,:) = ( r_CO2_uc(i,1,:)*mC + r_CO2_uc(i,2,:)*mOx + r_CO2_uc(i,3,:)*mOx )/(mOx + mOx + mC)
            vcm_CO2(i,:) = ( CO2_vxyz(i,1,:)*mC + CO2_vxyz(i,2,:)*mOx + CO2_vxyz(i,3,:)*mOx )/(mOx + mOx + mC)
         end do

         if (n_h2o > 0) then
         call msd_samp(msd_h2o,rcm_h2o_uc)
         call vcf_samp(vcf_h2o,vcm_h2o)
         end if
         if (n_gas > 0) then
         call msd_samp(msd_gas,rcm_GAS_uc)
         call vcf_samp(vcf_gas,vcm_GAS)
         end if
         call msd_samp(msd_O2,rcm_O2_uc)
         call vcf_samp(vcf_O2,vcm_O2)
         call msd_samp(msd_N2,rcm_N2_uc)
         call vcf_samp(vcf_N2,vcm_N2)
         call msd_samp(msd_CO2,rcm_CO2_uc)
         call vcf_samp(vcf_CO2,vcm_CO2)
      end if
      end if



      if (mod(istep,nsamp) == 0) then
      if (istep > nequib) then
         ! Kinetic energy of each species
         Ek = kinetic_energy(natom,vxyz,mass)
         Ek_h2o = kinetic_energy(n_h2o,n_atom_h2o,h2o_vxyz,h2o_mass)
         Ek_gas = kinetic_energy(n_GAS,n_atom_gas,GAS_vxyz,GAS_mass)
         Ek_O2 = kinetic_energy(n_O2,2,Ox_vxyz,Ox_mass)
         Ek_N2 = kinetic_energy(n_N2,2,N2_vxyz,N2_mass)
         Ek_CO2 = kinetic_energy(n_CO2,3,CO2_vxyz,CO2_mass)

         Ek_tot = Ek + Ek_O2+ Ek_N2 + Ek_CO2 + Ek_gas + Ek_h2o

         Ndf = 3*(natom + n_atom_gas*n_GAS + n_atom_h2o*n_h2o + 2*n_O2 + 2*n_N2 + 3*n_CO2) - 3

         T_kinetic = 2.0_wp*EK_tot/(K_ev*Ndf)

         write(io_energy,'(8g18.8)') time, Ek_tot, energy + ukeating + repulenergy*0.5_wp, &
            Ek_tot + energy + ukeating + repulenergy*0.5_wp
!         write(io_temp,'(9g18.8)') time, t_kinetic, &
!            2.0_wp*ek/(k_ev*(3*natom)), &
!            2.0_wp*Ek_h2o/(k_ev*(3*(n_atom_h2o*n_h2o))), &
!            2.0_wp*Ek_gas/(k_ev*(3*(n_atom_gas*n_gas))), &
!            2.0_wp*ek_o2/(k_ev*(3*(n_O2*2))), &
!            2.0_wp*ek_n2/(k_ev*(3*(n_N2*2))), &
!            2.0_wp*ek_co2/(k_ev*(3*(n_CO2*3)))

         call flush(io_energy);call flush(io_temp)
      end if
      end if

      call cpu_time(tp(4)); sumtp(3) = sumtp(3) + tp(4) - tp(3)
      if (mod(istep,nsamp) == 0) then
         write(*,'(i9,8g18.8)') istep,time,tp(2) - tp(1),tp(3) - tp(2),tp(4) - tp(3), tp(4) - tp(1), tp(4) - tp(0)
      end if

      if (make_movie) then
      if (mod(istep,nfreq_movie) == 0) then
         iframe = iframe + 1
         call new_frame_file(imve,'frame',iframe)
         call write_box_frame(imve,n_CO2,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
         call new_frame_file(imve,'frame_vel',iframe)
         call write_vel(imve,n_CO2,print_all = (mod(istep,nprint_all) == 0))
         call flush(imve)
      end if
      end if

!===================================================== end main loop
      end do main_loop

      if (n_h2o > 0) then
      call  msd_print(msd_h2o,dtreal)
      call  vcf_print(vcf_h2o,dtreal)
      end if
      if (n_gas > 0) then
      call  msd_print(msd_gas,dtreal)
      call  vcf_print(vcf_gas,dtreal)
      end if
      call  msd_print(msd_o2,dtreal)
      call  vcf_print(vcf_o2,dtreal)
      call  msd_print(msd_n2,dtreal)
      call  vcf_print(vcf_n2,dtreal)
      call  msd_print(msd_co2,dtreal)
      call  vcf_print(vcf_co2,dtreal)

      do i = 1,4
         print *,i,sumtp(i)
      end do
      close(io_energy,status='keep')
      close(io_temp,status='keep')

      call get_free_file_unit(io_config_final)
      open(unit=io_config_final,file='final.xyz')
      call write_xyz(io_config_final)
      close(io_config_final,status='KEEP')

CONTAINS

   FUNCTION frac_dens_wt()
      real(wp):: frac_dens_wt
      real(wp):: mass
      integer:: i
      mass = 0.0_wp
      do i = 1,natom
         select case(atom(i))
         case(iSilicon)
            mass = mass + mSi
         case(iOxygen)
            mass = mass + mOx
         case(iOxygenH)
            mass = mass + mOH
         case default
            stop 'unknown atom type'
         end select
      end do
      frac_dens_wt = mass/(8000*mSi + 16000*mOx)
   END FUNCTION frac_dens_wt

END PROGRAM simbox

