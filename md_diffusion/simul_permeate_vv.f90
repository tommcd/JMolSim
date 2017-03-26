
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
include 'frames_mod.f90'
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

PROGRAM SIMBOX
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
      USE frames_mod
      USE comvel_mod
      USE force_keating_mod
      USE vverlet_mod
      USE keating_parameters
      USE nlist_mod
      USE repul_energy_mod
      USE charges_mod
      USE insert_mod
      implicit none
      integer,parameter:: natom_unit=24
      integer,parameter:: nunitc = 4
      integer,parameter:: number_moves = 100000
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: mOH = 17.00_wp
      real(wp),parameter:: mN = 14.00_wp
      real(wp),parameter:: mC = 12.00_wp
      real(wp),parameter:: dt_fs = 0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,dtreal,qj,sumqj,rxyz12(3)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp):: CL,bl1,bl2,r3l,time,Ek,rad_Ox(3),rad_N2(3),rad_CO2(3),ULJEL
      real(wp),allocatable:: massi(:,:),Ox_massi(:,:,:),N2_massi(:,:,:),CO2_massi(:,:,:)
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic = 0,natomc,ix,iy,iz,istep
      integer:: itemp,itmp
      real(wp),allocatable:: coprxyz(:,:)
      character(32):: ctmp,c6*6,c5*5
!
      open(unit=20,file='connect_struct.out')
      open(unit=14,file='simbox.pdb')
      open(unit=15,file='neighbor.out')
      open(unit=17,file='energy.out')
!
      natom_max = 19900
!      CL = 0.7168_wp*0.162_wp/0.1552_wp
!      boxl = nunitc*CL
      boxl = 7.13286_wp
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      rad_Ox = (/sig_O_O2,sig_O_O2,0.0_wp/)*0.5_wp
      rad_N2 = (/sig_N_N2,sig_N_N2,0.0_wp/)*0.5_wp
      rad_CO2 = (/sig_C_CO2,sig_O_CO2,sig_O_CO2/)*0.5_wp
      n_O2 = 0!40
      n_N2 = 0!40
      n_CO2 = 0!40
      natom = natom_max
      T(1) = 15000.0_wp
      write (*,*) natom_max
!
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(Ox_xyz(1:n_O2,1:3,1:3),N2_xyz(1:n_N2,1:3,1:3),CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(Imige_Ox_xyz(1:n_O2,1:3,1:3),Imige_N2_xyz(1:n_N2,1:3,1:3),Imige_CO2_xyz(1:n_CO2,1:3,1:3))
      allocate(fxyz(1:natom_max,3),Ox_fxyz(1:n_O2,1:3,1:3),N2_fxyz(1:n_N2,1:3,1:3),CO2_fxyz(1:n_CO2,1:3,1:3))
      allocate(vxyz(1:natom_max,3),Ox_vxyz(1:n_O2,1:3,1:3),N2_vxyz(1:n_N2,1:3,1:3),CO2_vxyz(1:n_CO2,1:3,1:3))
      allocate(atom(1:natom_max),Ox_atom(1:n_O2,3),N2_atom(1:n_N2,3),CO2_atom(1:n_CO2,3))
      allocate(charge(natom),Ox_gas_charge(n_O2,3),N2_gas_charge(n_N2,3),CO2_gas_charge(n_CO2,3))
      allocate(proximity(natom_max,4))
      allocate(massi(natom_max,3),Ox_massi(1:n_O2,3,3),N2_massi(1:n_N2,3,3),CO2_massi(1:n_CO2,3,3))
      allocate(rxyzc(1:natom_unit,3),mol_vec(3,3))
      proximity = 0

      do i = 1,natom
         read (14,'(a6)',advance = 'no') c6
         read (14,*) itmp,ctmp,itmp,(rxyz(i,j),j = 1,3)
print '(i,a,i,3f12.6)', i,ctmp,itmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
print *,i,atom(i),'   ',atom_name(atom(i))
         if (atom(i) == iSilicon) massi(i,:) = 1.0_wp/mSi
         if (atom(i) == iOxygen)  massi(i,:) = 1.0_wp/mOx
         if (atom(i) == iOxygenH)  massi(i,:) = 1.0_wp/mOH
      end do

      rxyz(1:natom,:) = rxyz(1:natom,:)/angstrom
proximity = 0
      do i = 1,natom
         read(14,'(a32)') ctmp
print *,i,trim(ctmp)
         c5 = ctmp(6 + 1:6 + 5)
         read( unit=c5,fmt=* ) ii
         do j = 1,ncmax(atom(ii))
            c5 = ctmp(6 + 5*(j) + 1:6 + 5*(j) + 5)
            read( unit=c5,fmt=* ) proximity(ii,j)
print *,j,proximity(ii,j)
         end do
      end do

      call assign_charge(1,natom)

write(*,*) 'sum charge = ',sum(charge(1:natom))
!write(*,*) 'sum charge/qo = ',sum(charge(1:natom))/qi(iOxygen)
!write(*,*) 'n Si = ',count(atom==iSilicon)
!write(*,*) 'n O  = ',count(atom==iOxygen)
!write(*,*) 'n OH = ',count(atom==iOxygenH)
!write(*,*) 'ntot = ',count(atom==iSilicon) + count(atom==iOxygen) + count(atom==iOxygenH)
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



!      call INIT_NLIST(boxl,boxl,boxl,5.0_wp/angstrom)
!      call NEW_NLIST


open(66,file='oxygen.xyz')
write(66,*) natom + 2*n_o2 +2*n_N2 +3*n_CO2 + 2
write(66,*)'Silica + Gases'
do i = 1,natom
write(66,'(a,3f12.4)') atom_name(atom(i)),rxyz(i,:)*10.0_wp
end do

! ================== insertion procces==================================
!      mol_vec(1,1) = 0.05845_wp
!      mol_vec(1,2) = -0.05845_wp
!      mol_vec(1,3) = 0.0_wp
!      mol_vec(2,:) = 0.0_wp
!      mol_vec(3,:) = 0.0_wp
!   call insert_mol_N2_O2(n_O2,Ox_xyz,1000000,rad_Ox)
!
!do i = 1,n_O2
!write(66,'(a,3f12.4)')'O2 ',Ox_xyz(i,1,:)*10.0_wp
!write(66,'(a,3f12.4)')'O2 ',Ox_xyz(i,2,:)*10.0_wp
!end do
!
!      mol_vec(1,1) =-0.05485_wp
!      mol_vec(1,2) = 0.05485_wp
!      mol_vec(1,3) = 0.0_wp
!      mol_vec(2,:) = 0.0_wp
!      mol_vec(3,:) = 0.0_wp
!   call insert_mol_N2_O2(n_N2,N2_xyz,1000000,rad_N2)
!
!
!do i = 1,n_N2
!write(66,'(a,3f12.4)')'N2 ',N2_xyz(i,1,:)*10.0_wp
!write(66,'(a,3f12.4)')'N2 ',N2_xyz(i,2,:)*10.0_wp
!end do
!
!
!      mol_vec(1,1) = 0.0_wp
!      mol_vec(1,2) = 0.1162_wp
!      mol_vec(1,3) = -0.1162_wp
!      mol_vec(2,:) = 0.0_wp
!      mol_vec(3,:) = 0.0_wp
!   call insert_mol_CO2(n_CO2,CO2_xyz,1000000,rad_CO2)
!
!do i = 1,n_CO2
!write(66,'(a,3f12.4)')'CO ',CO2_xyz(i,1,:)*10.0_wp
!write(66,'(a,3f12.4)')'OC ',CO2_xyz(i,2,:)*10.0_wp
!write(66,'(a,3f12.4)')'OC ',CO2_xyz(i,3,:)*10.0_wp
!end do
!open(66,file='oxygen.xyz')
!write(66,*) natom + 2*n_o2 +2*n_N2 +3*n_CO2 + 2
!write(66,*)'Silica + Gases'
!do i = 1,natom
!write(66,'(a,3f12.4)') atom_name(atom(i)),rxyz(i,:)*10.0_wp
!end do
!write(66,'(a,3f12.4)')'Ar ',(/boxl2,boxl2,boxl2/)*10.0_wp
!write(66,'(a,3f12.4)')'Ar ',-(/boxl2,boxl2,boxl2/)*10.0_wp

!!========================================================================
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


      close(15,status='keep')
      Ox_vxyz = 0.0_wp
      N2_vxyz = 0.0_wp
      CO2_vxyz = 0.0_wp
      call COMVEL(natom,n_CO2,n_O2,n_N2,T(1),vxyz,Ox_vxyz,N2_vxyz,CO2_vxyz)

      vxyz = vxyz*sqrt(massi)
      Ox_vxyz = Ox_vxyz*sqrt(Ox_massi)
      N2_vxyz = N2_vxyz*sqrt(N2_massi)
      CO2_vxyz = CO2_vxyz*sqrt(CO2_massi)
      call set_bond_angle_lists()

 do i = 1,natom
   vxyz(i,2:3) = 0.0
   vxyz(i,1) = vxyz(12,1)
 end do

 rxyz12 = rxyz(12,:)

      fxyz = 0.0_wp
      Ox_fxyz = 0.0_wp
      N2_fxyz = 0.0_wp
      CO2_fxyz = 0.0_wp
      energy = 0.0_wp
      call force_keating
!      call repul_energy(1,natom)
!      call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!      call O2_stretching
!      call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!      call N2_stretching
!      call FORCE_ENERGY_CO2_LJ_EL(1,natom,1,n_CO2,ULJEL)
!      call CO2_stretching
!      call CO2_angle_bending

      time = 0.0_wp
      istep = 0
      dtreal = 1.0_wp
      dt = dtreal*dt_fs
         do n = 1, number_moves
         write(*,*) n
            istep = istep + 1
               coprxyz = rxyz

            call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
!            forall(i=1:natom) rxyz(i,1:3) =rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)
              rxyz(12,:) = rxyz12
!            call gas_vverlet_a(dt,Ox_massi,Ox_xyz,Ox_vxyz,Ox_fxyz)
!            forall(i=1:n_O2,j=1:3) Ox_xyz(i,j,1:3) =Ox_xyz(i,j,1:3)-boxl*anint(Ox_xyz(i,j,1:3)*boxli)
!
!            call gas_vverlet_a(dt,N2_massi,N2_xyz,N2_vxyz,N2_fxyz)
!            forall(i=1:n_N2,j=1:3) N2_xyz(i,j,1:3) =N2_xyz(i,j,1:3)-boxl*anint(N2_xyz(i,j,1:3)*boxli)
!
!            call gas_vverlet_a(dt,CO2_massi,CO2_xyz,CO2_vxyz,CO2_fxyz)
!            forall(i=1:n_CO2,j=1:3) CO2_xyz(i,j,1:3) =CO2_xyz(i,j,1:3)-boxl*anint(CO2_xyz(i,j,1:3)*boxli)

            ! check if the neighbour list needs to be updated
!                  do i = 1,natom
!               if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3))) then
!                  call NEW_NLIST
!                  EXIT
!               end if
!                  end do

            fxyz = 0.0_wp
            Ox_fxyz = 0.0_wp
            N2_fxyz = 0.0_wp
            CO2_fxyz = 0.0_wp
            energy = 0.0_wp
            call force_keating
!            call repul_energy(1,natom)
!            call FORCE_ENERGY_LJ_EL_O2(1,natom,1,n_O2,ULJEL)
!            call O2_stretching
!            call FORCE_ENERGY_LJ_EL_N2(1,natom,1,n_N2,ULJEL)
!            call N2_stretching
!            call FORCE_ENERGY_CO2_LJ_EL(1,natom,1,n_CO2,ULJEL)
!            call CO2_stretching
!            call CO2_angle_bending

            call vverlet_b(dt,massi,vxyz,fxyz)
!            call gas_vverlet_b(dt,Ox_massi,Ox_vxyz,Ox_fxyz)
!            call gas_vverlet_b(dt,N2_massi,N2_vxyz,N2_fxyz)
!            call gas_vverlet_b(dt,CO2_massi,CO2_vxyz,CO2_fxyz)
              time = time + dtreal
            if (mod(istep,1) == 0) then
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
!              do m=1, n_N2
!                do i=1,2
!              Ek = Ek + (0.5_wp/N2_massi(m,i,1))*dot_product(N2_vxyz(m,i,:),N2_vxyz(m,i,:))
!                end do
!              end do
!              do m=1, n_CO2
!                do i=1,3
!              Ek = Ek + (0.5_wp/CO2_massi(m,i,1))*dot_product(CO2_vxyz(m,i,:),CO2_vxyz(m,i,:))
!                end do
!              end do


             write (17,*) time, Ek + energy ! + repulenergy/2.0_wp
!             write (999,'(8f12.6)') Ek+energy +repulenergy/2.0_wp,Ox_xyz(1,1,:),Ox_xyz(1,2,:)
            call flush(17)
            end if
         end do
      close(17,status='keep')

      !call write_xyz(20,0)
      call write_frame(20,1,natom)
      close(14,status='KEEP')



end program simbox


