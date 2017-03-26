
include 'precision_mod.f90'
include 'math_const_mod.f90'
include 'global_vars.f90'
include 'sort_mod.f90'
include 'files_mod.f90'
include 'seaton_mod.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'strip_cell_mod.f90'
!include 'charges_mod.f90'
include 'bond_list_mod.f90'
include 'frames_mod.f90'
include 'rand_mod.f90'
include 'matvec3d_mod.f90'
include 'nlist_mod.f90'
include 'verlet_list_nl.f90'
include 'simbox_mod2.f90'

!include 'tetra_coords_mod.f90'
!include 'quat2mat_mod.f90'
!include 'ran_point_sphere_mod.f90'
!include 'Keating_parameters_Burlakov.f90'
!include 'repul_energy_VL_Burlakov.f90'
include 'Keating_parameters_vonAlfthan.f90'
include 'repul_energy_VL_vonAlfthan.f90'
include 'Keating.f90'
include 'HKNonLattice.f90'
!include 'probe_mol_mod.f90'
!include 'voidage4_mod.f90'
!include 'deposition4_mod_H.f90'
include 'rdf_mod.f90'
include 'bin_mod.f90'
include 'bond_angle_mod.f90'
include 'relax2_VL_mod_2.f90'
!include 'select_atom_mod.f90'
!include 'attachment_NL_mod_H.f90'
!include 'add_water2_VL_mod_H.f90'
!include 'relax_dangling2_VL_NL_mod_H.f90'
!include 'surface_atoms_mod.f90'
!include 'lj_el_mod.f90'
!include 'Henrys_law_calc_mod.f90'


PROGRAM unit_cell_relax
      USE precision_mod, only: wp
      USE global_vars
      USE constants_mod
      USE atom_types
      USE connectivity_mod
      USE coordinates_mod
      USE rand_mod
      USE frames_mod
      USE simbox_mod
      USE Relax_mod
      USE nlist_mod
      USE verlet_list
      USE HKNonLattice_mod
      USE rdf_mod
      USE Keating_parameters, only: ASiO
      USE files_mod
      USE Relax_mod
      USE keating,only:energy4
      USE repul_energy_mod
      USE seaton_mod
      USE bond_list_mod
      USE bond_angle_mod
      USE bin_mod
      implicit none
      logical,parameter:: make_movie = .TRUE.
      real(wp),parameter:: aSiOSi = pi*(140.0_wp/180.0_wp)
      real(wp),parameter:: r_ox = 0.15_wp
      real(wp),parameter:: d_ox = 2.0_wp*r_ox
      real(wp):: rverlet!,top_atom,void_crit,dlayer !!,z_pos
      integer:: i,ii,jj,kk,mm,L,imve,iu!,N_OH_groups,ntotal,iz,ipconfig
      integer:: NR_TOT,NR,NPRINT,iu1,iu2,iu3,ntypes,na
      real(wp):: timep(0:10),delt,delt0,Energy!,t0,t1
      character:: date*8,time*10
      logical:: quench,scale_repul,fexists,fopen!,success,ok
      !integer:: ntrial_henry
      !real(wp):: khenry,facc,fover,uljel

      OPEN(unit = 10,FILE = 'unit_cell.dat',STATUS = 'OLD')
!-----------------------------------------------------------
!     Read in simulation data
!     nattached = number of molecules already attached
!     ntotal = total number of molecules to be attached
!     nseed = number of initial surface oxygens
!     boxl = box side length of cell
!     ipconfig = frequency of xyz,RDF dump
!     irand_seed = random number seed
!-----------------------------------------------------------

      READ(10,*) Etot
      READ(10,*) delt0
      READ(10,*) irand_seed
      READ(10,*) rverlet
      READ(10,*) NR_TOT
      READ(10,*) NPRINT
      close(UNIT = 10,STATUS = 'KEEP')
      
      quench = .TRUE.
      bin_width = 0.2_wp
      max_dev = 71.0_wp
      ntypes = 2
      na = 0

      
      CALL simbox(rverlet)
      CALL NEW_VLIST 
      proximity = 0
      NR = NR_TOT/NPRINT
      nattached = 0
      delt = delt0

      OPEN(unit= 15,file= 'cell.xyz')
      OPEN(unit= 19,file= 'prox_data.out')

      CALL cpu_time(timep(0))
      CALL read_proximity(19,nattached)
      
      CALL init_bins(bin_width,max_dev,ntypes)
      CALL init_bond_angle(bin_width,max_dev,ntypes)
      
!      CALL new_bond_file(iu,'init_anal',0) 
!      CALL bond_angle(1,natom,ntypes)
!      do mm =  lbound(bin_med_angle,1), ubound(bin_med_angle,1) 
!      !write(iu,'(i4,2(2x,i8))') mm, abin(mm,0), abin(mm,1) 
!      write(iu,'(f9.5,2x,i8)') bin_med_angle(mm,0), abin(mm,0)
!      end do       
!      close(iu, STATUS = 'KEEP')


      CALL new_bond_file(iu1,'O_temp',0)
      CALL new_bond_file(iu2,'Si_temp',0)
      CALL new_bond_file(iu3,'all_temp',0)  

      print*, 'iu1 = ', iu1
      print*, 'iu2 = ', iu2
      print*, 'iu3 = ', iu3
!      CALL bond_angle(1,natom,ntypes)
!      do kk =  lbound(bin_med_angle,1), ubound(bin_med_angle,1) 
!      !write(iu,'(i4,2(2x,i8))') kk, abin(kk,0), abin(kk,1) 
!      write(iu,'(f9.5,2x,i8)') bin_med_angle(kk,0), abin(kk,0)
!      end do       
!      close(iu, STATUS = 'KEEP')
      
! the following should produce a relaxation at max temp for first 
! NR steps allowing it to melt and is then quenched reducing the 
! temperature
      do i=0,NPRINT+1
         if (i<=1) then
            quench = .FALSE.
            Energy = delt
         else 
            quench = .TRUE.
!            Energy = (NPRINT-i)*Etot

!            ! for a linear cooling rate use the following
!            Energy = ((real(NPRINT-(i-1)))/real(NPRINT))*delt0
!            delt = ((real(NPRINT-(i-2)))/real(NPRINT))*delt0 - (1/NR)*delt0
            ! for an essentially quadratic rate then use
            delt = Energy - (1/NR)*delt0
            Energy = 0.8*Energy
            
         end if 
   
         ! the Etot and delt are the lower and upper (resp.) bounds for the 
         ! relazing temperatures.
         CALL relax(NR,1,natom,Energy,delt,quench)
         CALL new_frame_file(imve,'frame',i)
         CALL write_frame(imve,1,natom)

         ! after each relaxation step the bond angles are calculated and written
         ! to separate files for each respective angle type
         ! Si-O-Si bonds
         CALL bond_angle(1,natom,ntypes)
         
         inquire(iu1,exist= fexists,opened = fopen)
         if (fopen) print*, 'ok'

         
         !CALL new_bond_file(iu,'O_temp',i+1)
         do ii = lbound(bin_med_angle,1), ubound(bin_med_angle,1)
!   print *,'ii = ',ii
!   print *,'bin_med_angle = ', bin_med_angle
!   print *, lbound(bin_med_angle,1), ubound(bin_med_angle,1)
!   print *, lbound(bin_med_angle,2), ubound(bin_med_angle,2)
!   write(iu1,'(f9.5)') bin_med_angle(ii,0)
!   write(iu1,'(f9.5)') bin_med_angle(ii,1)
!   write(iu1,'(i8)') abin(ii,0)
         write(iu1,'(f9.5,2x,i8)') bin_med_angle(ii,0), abin(ii,0)
         end do
         write(iu1,'(/)')
                         
         ! O-Si-O bonds 
         !CALL new_bond_file(iu,'Si_temp',i+1)
         do jj =  lbound(bin_med_angle,1), ubound(bin_med_angle,1) 
         write(iu2,'(f9.5,2x,i8)') bin_med_angle(jj,1), abin(jj,1) 
         end do
         write(iu2,'(/)')         
     
! alternatively the angles could be both written to the one file and have the 
! x-value as the deviation rather than actual value of angles
         !CALL new_bond_file(iu,'temp',i+1)
         do kk =  lbound(bin_med_angle,1), ubound(bin_med_angle,1) 
         write(iu3,'(i4,2(2x,i8))') kk, abin(kk,0), abin(kk,1) 
         end do 
         write(iu3,'(/)')         

         
!      do i=1, NPRINT
!print*, 'i (nprint loop) = ', i
!! the Etot and delt are the lower and upper (resp.) bounds for the 
!! relazing temperatures.
!         call relax(NR,1,natom,Etot,delt,alist,na,quench,scale_repul)
!         call new_frame_file(imve,'frame',i)
!         call write_frame(imve,1,natom)
!      end do 

      end do 
      call cpu_time(timep(4))
      call date_and_time(date,time)
      print*, ''
      print *, 'Simulation time: ', timep(4)-timep(0)
      
      close(iu1,STATUS = 'KEEP') 
      close(iu2,STATUS = 'KEEP')  
      close(iu3,STATUS = 'KEEP') 
      close(15)
      close(19)

      OPEN(unit= imve)
      CALL pdb2xyz_mod(imve)
      CLOSE(imve)

      ! initial translation of box occurs as follows 
      rxyz(1:natom,1:2) = rxyz(1:natom,1:2) + boxl2

      ! image in y direction
      rxyz(natom+1:2*natom,2) = rxyz(1:natom,2) - boxl
      atom(natom+1:2*natom) = atom(1:natom)      
      
      ! image in x direction 
      rxyz(2*natom+1:4*natom,1) = rxyz(1:2*natom,1) - boxl
      atom(2*natom+1:4*natom) = atom(1:2*natom)
      
      ! image in z direction
      rxyz(4*natom+1:8*natom,3) = rxyz(1:4*natom,3) + boxl
      rxyz(8*natom+1:12*natom,3) = rxyz(1:4*natom,3) - boxl
      atom(4*natom+1:8*natom) = atom(1:4*natom)     
      atom(8*natom+1:12*natom) = atom(1:4*natom)




      

END PROGRAM  unit_cell_relax

