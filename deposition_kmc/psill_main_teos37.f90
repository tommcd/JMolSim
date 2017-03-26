
INCLUDE 'precision_mod.f90'
INCLUDE 'rand_mod.f90'
INCLUDE 'global_vars_mod_KMC.f90'
INCLUDE 'M_intcat_jlg.f90'
INCLUDE 'misc_legacy_mod.f90'
INCLUDE 'math_const_mod.f90'
INCLUDE 'phys_const_mod.f90'
INCLUDE 'command_line_mod.f90'
INCLUDE 'command_arg_mod.f90'
INCLUDE 'files_mod.f90'
INCLUDE 'sort_mod.f90'
INCLUDE 'global_vars_mod_KMC.f90'

INCLUDE 'list_mod.f90'
INCLUDE 'atom_types_mod_teos.f90'
INCLUDE 'coordinates_mod.f90'
INCLUDE 'nlist_mod_zplus.f90'
INCLUDE 'connectivity_teos_mod.f90'
INCLUDE 'check_structure_mod_TEOS.f90'
INCLUDE 'math_const_mod.f90'
INCLUDE 'phys_const_mod.f90'
INCLUDE 'Dreiding_parameters_mod.f90'
INCLUDE 'seaton_mod.f90'
INCLUDE 'charmm_mod.f90'
INCLUDE 'Keating_parameters_mod_vonAlfthan.f90'
INCLUDE 'constants_mod.f90'
INCLUDE 'teos_atomlst_mod.f90'
INCLUDE 'find_atoms_mod.f90'
INCLUDE 'charges2_mod_teos.f90'
INCLUDE 'bond_angle_types_mod_teos.f90'
INCLUDE 'bond_angle_list_mod.f90'
INCLUDE 'box_frames_mod.f90'
INCLUDE 'frames_mod.f90'
INCLUDE 'verlet_list_nl.f90'
INCLUDE 'force_keating_mod.f90'
INCLUDE 'repul_force_NL_vonAlfthan.f90'
INCLUDE 'energy_keating_mod.f90'
INCLUDE 'force_en_el_mod.f90'
INCLUDE 'lj_el_mod_O2_N2_CO2.f90'
INCLUDE 'force_energy_minimize_mod.f90'
INCLUDE 'relax_mc_mod.f90'
INCLUDE 'hydrogens_mod.f90'
INCLUDE 'bond_list_mod.f90'
INCLUDE 'matvec3d_mod.f90'
INCLUDE 'rotate_axis_mod.f90'
INCLUDE 'insert_seeds_mod_H.f90'
INCLUDE 'tetra_coords_mod.f90'
INCLUDE 'quat2mat_mod.f90'
INCLUDE 'ran_point_sphere_mod.f90'
INCLUDE 'HKNonLattice2.f90'
INCLUDE 'datastore_mod.f90'
INCLUDE 'probe_mol_mod_teos.f90'
INCLUDE 'lj_el_mod.f90'
INCLUDE 'Henrys_law_calc_mod.f90'
INCLUDE 'voidage4_mod.f90'
INCLUDE 'deposition4_mod_H.f90'
INCLUDE 'rdf_mod.f90'
INCLUDE 'relax_VL_mod.f90'
INCLUDE 'attachment_mod_teos.f90'
INCLUDE 'add_hydrogen_mod.f90'
INCLUDE 'glassbond_mod.f90'
INCLUDE 'dissociation_mod.f90'
INCLUDE 'water_attack_VL_mod_H.f90'
INCLUDE 'condensation_mod_H.f90'
INCLUDE 'surface_atoms_mod.f90'
INCLUDE 'switch.f90'
INCLUDE 'relax_bonds_VL_NL_mod.f90'
INCLUDE 'list_atoms_in_cells_mod.f90'
INCLUDE 'reaction_rate_mod.f90'


PROGRAM main_deposition
      USE precision_mod, only: wp
      USE global_vars_mod
      USE atom_types_mod
      USE rand_mod
      USE teos_atomlst_mod
      USE math_const_mod
      USE phys_const_mod
      USE coordinates_mod
      USE connectivity_mod
      USE check_structure_mod
      USE seaton_mod
      USE Dreiding_parameters_mod
      USE charmm_mod
      USE Keating_parameters_mod
      USE constants_mod
      USE bond_angle_types_mod
      USE bond_angle_list_mod
!     USE BOX_FRAMES_MOD
      USE FRAMES_MOD
      USE charges_mod
      USE nlist_mod
      USE verlet_list_mod
      USE force_keating_mod
      USE repul_energy_mod
      USE force_en_el_mod
      USE energy_keating_mod
      USE lj_el_old_mod
      USE lj_el_mod
      USE FORCE_ENERGY_minimize_mod
      USE relax_mc_mod
      USE list_mod
      USE bond_list_mod
      USE hydrogens_mod
      USE rotate_axis_mod
      USE HKNonLattice_mod
      USE tetra_coords_mod
      USE add_hydrogen_mod
      USE ran_point_sphere
      USE quat2mat_mod
      USE find_atoms_mod
      USE probe_mol_mod
      USE voidage_mod
      USE deposition_mod
      USE glassbond_mod
      USE rdf_mod
      USE water_attack_mod
      USE files_mod
      USE relax_mod
      USE attachment_mod
      USE insert_seeds_mod
      USE surface_atoms_mod
      USE Henrys_law_calc_mod
      USE condensation_mod
      USE datastore_mod
      USE dissociation_reaction_mod
      USE relax_bonds_mod
      USE list_atoms_in_cells_mod
      USE reaction_rate_mod

      implicit none
      real(wp),parameter:: ulength = 1.0e-10_wp*angstrom
      real(wp),parameter:: uenergy = 1.60217733e-19_wp
      real(wp),parameter:: upressure = uenergy/(ulength**3)
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: aSiOSi = pi*(140.0_wp/180.0_wp)
      real(wp):: top_atom,void_crit,rverlet
      integer:: i,j,ipconfig,ntotal,imve,ib,ir,k,ier
      real:: timep(0:10)
      character:: date*8,ctime*10
      logical:: success,ok
!
! Kinetic Monte Carlo

real(wp):: psum
type(list):: alst(16),oh_lst,ohb_lst
real(wp):: time,delta_time,xi1
real(wp),parameter:: torr2Pa = 101325.0_wp/760.0_wp
integer:: nr,L,iSi,iu
real(wp):: Temp_K,kF(8) = 0.0_wp,A(8) = 0.0_wp,Ea(8) = 0.0_wp
real(wp),parameter:: R_gas_cal = 1.987_wp
!
      call data_file
 !
      print '(a,g18.9,a)','ulength   = ',ulength,' m'
      print '(a,g18.9,a)','uenergy   = ',uenergy,' J'
      print '(a,g20.10,a)','upressure = ',upressure,' Pa'
      print '(a,g18.9,a,g18.9,a)','P_SiOC2H5 = ',P_SiOC2H54,' = ',P_SiOC2H54*upressure,' Pa'
      print '(a,g18.9,a)',        '        = ',P_SiOC2H54*upressure/torr2Pa,' torr'
      print '(a,g18.9,a)','kT   = ',kboltzT,' eV'
      print *,'T = ',Temp_K,' K'
      print '(a,g18.9)','exp(-EA/kT) = ',exp_Ea_kt

print *,'1st 3 random numbers, ',rand(),rand(),rand()
!
      natom_max = 2*nseed + 25*ntotal
      allocate(rxyz(1 - nseed:natom_max,3))
      allocate(atom(1 - nseed:natom_max))
      allocate(charge(1 - nseed:natom_max))
      allocate( proximity(natom_max,4) )
      allocate(copproximity(natom_max,size(proximity,dim=2)),copatom(natom_max))
      allocate(coprxyz(natom_max,3))
      call Init_Verlet_List(rverlet,0.26_wp)
      call Init_HKNonLattice(natom_max)
      call LJ_init_old
!
      proximity = 0
      boxli = 1.0_wp/boxl
      boxl2 = boxl/2.0_wp
!     boxl2n = -boxl2
!     boxl2i = 1.0_wp/boxl2
!
      open(unit=15,file='psil_config.xyz')
      open(unit=20,file='psil_topology.out')
      open(unit=14,file='psil_RDF.out',POSITION = 'APPEND')
      open(unit=21,file='psil_rand.out')
      open(unit=22,file='psil_dens.out')
     !open(unit=23,file='psil_void_SiO4.out')
      open(unit=26,file='psil_void_H2O.out')
      open(unit=24,file='psil_topatom.out')
      open(unit=25,file='psil_time.out')
      open(unit=29,file='psil_kmc_history.out')
      open(unit=77,file='psil_SiO4_attempt.out')
      open(unit=78,file='psil_SiO4_add.out')
      open(unit=79,file='psil_H2O_add.out')
      open(unit=80,file='psil_H2O_remove.out')
      open(unit=81,file='psil_bond_switch.out')

      if (kmc_step == 0) then
         write(*,*) 'starting from scratch'
         call insert_seeds
         n_cluster = nseed
         n_cluster_old = n_cluster
         forall(i = 1:n_cluster) atomL(i) = i  ! Label each seed as a cluster
         if (make_movie) then ! write out first frame of movie
            call new_frame_file(imve,'frame',0)
            call write_frame(imve,1 - nseed,natom)
         end if
         top_atom = 0.0_wp
      else
         call read_xyz1(15,kmc_step)
         call read_proximity(20,kmc_step)
         call read_rand_state(21)
         call HKNonLattice(natom,proximity,n_cluster,atomL)
         n_cluster_old = n_cluster
         write(*,*) 'n_cluster = ',n_cluster
         write(*,*) 'atomL = '; write(*,'(10i6)') atomL(1:natom)
         top_atom = maxval(rxyz(1:natom,3))
         forall(i = 1:nseed) rxyz(i - nseed,:) = rxyz(i,:) - (/0.0_wp,0.0_wp,bondl_SiO/)
         atom(1 - nseed:0) = iSilicon
      end if
!
      call init_probe_mols()
      call INIT_NLIST(boxl(1),boxl(2),500.0_wp/angstrom,max(5.0_wp/angstrom,2.0_wp*maxval(probe_tip3p%rad)),natom_max)
      call NEW_NLIST(1,natom)
      call new_vlist
      call cpu_time(timep(0))
!
!----------------------------------------------------------------------
      time = 0.0_wp
      main_loop: do while (kmc_step < ntotal)
      call cpu_time(timep(1))
      kmc_step = kmc_step + 1
!
      call assign_charge(1,natom)
!print*,'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
      ! write(*,*)'net charge = ',sum(charge(1 - nseed:natom))
!
      dlayer = (5.0_wp/angstrom)
      volbin = boxl(1)*boxl(2)*dlayer
      call henry_profile(pr_TEOS,dlayer,1)
      call henry_profile(probe_tip3p,dlayer,2)
      call den_profile_rxn(1,natom,dlayer)
      call prob_dist(henryk(1:nbin,1),inthenryk(1:nbin,1),void_crit,itop_nonp_TEOS)
      call prob_dist(henryk(1:nbin,2),inthenryk(1:nbin,2),void_crit,itop_nonp_H2O)
print *,'reached - - -- - -- - -- - -- - -- - 1'
!print *,kF
print*,'nbin = ', nbin
100   CONTINUE
      call reaction_rate_cal(kF)
!
print *,'reached - - -- - -- - -- - -- - -- - 2'

if (rate_sum == 0.0_wp) STOP 'NO REACTION POSSIBLE'

      xi1 = rand()
      psum = 0.0_wp
      bin_loop: do ib = 1,nbin
      do ir = 1,nrxn
         psum = psum + rate(ir,ib)/rate_sum
         if (psum >= xi1) then
            exit bin_loop
         end if
      end do
      end do bin_loop
print *,'reached - - -- - -- - -- - -- - -- - 3'

!
!     carry out reaction #ir in bin #ib
      call check_proximity(1,natom,ok,ier)
            if (.NOT.ok) then
               write(*,*) 'proximity array is NOT consistent'
               write(*,*) 'check_proximity ',ok,ier
               STOP
            end if
            call check_proximity2_H(1,natom,atom(1:natom),proximity,ok,ier)
print*,'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
            if (.NOT.ok) then

              write(*,*) 'proximity array is NOT consistent'
              write(*,*) 'check_proximity2_H ',ok,ier

   do i = 1,natom
      k = count(proximity(i,:) > 0)
     write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
      i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
   end do
              STOP
            end if
         if (ANY(atom(proximity(1:nseed,1)) == iOxygen)) then
            write(*,*) 'seeds should not be bonded to Oxygens'
            STOP
         end if
         print *,'i am going inside the  randomly selected bin reaction bin'
      call reaction_in_bin (ir,ib,success)
      if (.NOT.success) then
         rate_sum = rate_sum - rate(ir,ib)
         rate(ir,ib) = 0.0_wp
         GOTO 100
      end if
      !write(29,'(g18.8,i4,2i9)') rate_sum,ir,ib,nbin
      write(29,*)rate_sum,ir,ib,nbin,henryk(ib,1),henryk(ib,2)
print *,'reached - - -- - -- - -- - -- - -- - 4'
!print *,'i am finished recation and one KMC step'
      call check_proximity(1,natom,ok,ier)
         if (.NOT.ok) then
            write(*,*) 'proximity array is NOT consistent'
            write(*,*) 'check_proximity',ok,ier
            STOP
         end if
      call check_proximity2_H(1,natom,atom(1:natom),proximity,ok,ier)
         if (.NOT.ok) then
           write(*,*) 'proximity array is NOT consistent'
           write(*,*) 'check_proximity2_H',ok,ier


!do i = 1,natom
!   k = count(proximity(i,:) > 0)
!  write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
!   i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
!end do
           STOP
         end if
         if (ANY(atom(proximity(1:nseed,1)) == iOxygen)) then
            write(*,*) 'seeds should not be bonded to Oxygens'
            STOP
         end if
!
!     advance the kmc clock
      delta_time = -log(rand())/rate_sum
      time = time + delta_time
!
      call cpu_time(timep(2))
      call new_vlist
!
      call RELAX_MC(NRELAX,1,natom, kboltzT)
!      call RELAX(NRELAX, nseed+1, natom, kboltzT)
!      print*,'999999999999999999999999999999999999999999999999999'
!      call adjust_hydrogens()
!
      call get_free_file_unit(iu)
      open(unit=iu, file = 'test1.xyz',status='UNKNOWN')
      write(iu,*) natom
      write(iu,'(a,i6)') 'test_file',natom
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)),(rxyz(i,j),j = 1,3)
      end do
      write(iu,*)
      close(iu)

      call cpu_time(timep(3))
!      call relax_bonds( count(atom(nseed+1:natom) == iOxygen) )
!      call adjust_hydrogens()
!
      call get_free_file_unit(iu)
      open(unit=iu, file = 'test2.xyz',status='UNKNOWN')
      write(iu,*) natom
      write(iu,'(a,i6)') 'test_file',natom
      do i = 1,natom
         write(iu,'(a2,3(1x,f14.8))') atom_name(atom(i)),(rxyz(i,j),j = 1,3)
      end do
      write(iu,*)
      close(iu)


      call cpu_time(timep(4))
      call date_and_time(date,ctime)
!
!     write(*,*) 'n O = ',count(atom(1:natom) == iOxygen)
!     write(*,*) 'n Si= ',count(atom(1:natom) == iSilicon)
      write(* ,'(i7,a10,2x,a10,6(3x,f0.3))') natom,date,ctime,timep(2) - timep(1), &
         timep(3) - timep(2),timep(4) - timep(3),timep(4) - timep(0)
      write(25,'(i7,a10,2x,a10,6(2x,f0.3))') natom,date,ctime,timep(2) - timep(1), &
         timep(3) - timep(2),timep(4) - timep(3),timep(4) - timep(0)
      top_atom = maxval(rxyz(1:natom,3))
      write(24,'(f12.6,2i7,f12.6)') time,kmc_step,natom,top_atom
!
      if (make_movie) then
         call new_frame_file(imve,'frame',kmc_step)
         call write_frame(imve,1,natom)
      end if
!
      if ( mod(kmc_step,ipconfig) == 0 .or. (kmc_step == ntotal) ) then
!--------Periodically check the connectivity and write out xyz coordinates,
!--------RDF, Connectivity, etc.
         call check_proximity(1,natom,ok,ier)
         if (.NOT.ok) then
            write(*,*) 'proximity array is NOT consistent'
            write(*,*) 'check_proximity ',ok,ier
            STOP
         end if
         call check_proximity2_H(1,natom,atom(1:natom),proximity,ok,ier)
         if (.NOT.ok) then
           write(*,*) 'proximity array is NOT consistent'
           write(*,*) 'check_proximity2_H ',ok,ier

!do i = 1,natom
!   k = count(proximity(i,:) > 0)
!  write(iu,'(i6,1x,a,1x,i5,4x,6(1x,i5))') &
!   i,atom_name(atom(i)),k,(proximity(i,j),j = 1,ncmax(atom(i)))
!end do
           STOP
         end if
         if (ANY(atom(proximity(1:nseed,1)) == iOxygen)) then
            write(*,*) 'seeds should not be bonded to Oxygens'
            STOP
         end if
         call rdf_calc
         call rdf_print(14,kmc_step)
         call write_xyz1(15,kmc_step)
         rewind(15)
         call write_proximity(20,kmc_step)
         rewind(20)
         call write_rand_state(21)
         rewind(21)
!         if (use_voidage_calc) call den_profile(1,natom,(5.0_wp/angstrom))
!         call write_den_profile(22)
!         call voidage_profile(probe_SiO4,dlayer)
!         call write_voidage(23)
!         call SURFACE_ATOMS(nseed+1,natom,(1.5_wp/angstrom))
!         write(*,*) 'fraction of surface atoms = ',fsurfa
!         !call write_surface_atoms_zbin(191)
!         call voidage_profile(probe_tip3p,      boxli = 1.0_wp/boxl

!         call write_voidage(26)
      end if

      end do main_loop
!----------------------------------------------------------------------
      close(14)
      close(15)
      close(20)
      close(21)
      close(22)
!     close(23)
      STOP
CONTAINS

     SUBROUTINE data_file
!!     real(wp):: kF(:),A(:),Ea(:)
     open(unit=10,file='kmc.dat',status='OLD')
!-----------------------------------------------------------
!     Read in simulation data
!     kmc_step = number of kmc steps so far
!     ntotal = total number of kmc steps
!     nseed = number of initial surface OH groups
!     boxl = box side length of cell
!     ipconfig = frequency of xyz,RDF dump
!     irand_seed = random number seed
!-----------------------------------------------------------
      READ(10,*) kmc_step  ;print *,'kmc_step = ',kmc_step
      READ(10,*) ntotal
      READ(10,*) nseed
      READ(10,*) Temp_K
      READ(10,*) P_SiOC2H54
      READ(10,*) fwater
      READ(10,*) boxl(1), boxl(2),boxl(3)
      READ(10,*) void_crit
      READ(10,*) ipconfig  ;print *,' ipconfig = ',ipconfig
      READ(10,*) irand_seed
      READ(10,*) NRELAX
      READ(10,*) rverlet  ;print *,' rverlet = ',rverlet

      close(unit=10,status='KEEP')
!
      open(unit=11,file='rateconst.dat',status='OLD')
      READ(11,*) exp_Ea_kt

      READ(11,*) nr
      do i = 1,nr
        READ(11,*) A(i),Ea(i)
        kF(i) = A(i)*exp(-Ea(i)/(R_gas_cal*Temp_K))
        print*,i, A(i),Ea(i),kF(i)
      end do
      close(unit=11,status='KEEP')
      kboltzT = K_ev*Temp_K
     ! exp_Ea_kt = exp(-del_rxn/kboltzT)
      P_SiOC2H54 = P_SiOC2H54/upressure
      P_H2O = fwater*P_SiOC2H54
      END SUBROUTINE data_file

END PROGRAM main_deposition

