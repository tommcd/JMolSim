
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
include 'bond_list_mod.f90'
include 'frames_mod.f90'
include 'comvel_mod.f90'
include 'force_keating.f90'
include 'Keating_parameters_vonAlfthan.f90'
include 'vverlet_mod.f90'
include 'nlist_mod.f90'
include 'repul_force_NL_vonAlfthan.f90'

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
      USE HKNonLattice_mod
      USE frames_mod
      use comvel_mod
      use force_keating_mod
      use vverlet_mod
      use keating_parameters
      use nlist_mod
      use repul_energy_mod
      implicit none
      integer,parameter:: natom_unit = 24
      integer,parameter:: nunitc = 4
      integer,parameter:: number_moves = 300000
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp),parameter:: dt_fs =0.982269e-2_wp
      real(wp):: r3(3),T(10),dt,dtreal
      integer,allocatable:: crossbond_x(:),crossbond_y(:),crossbond_z(:)
      real,allocatable:: rxyzc(:,:)
      integer,allocatable:: atomc(:)
      real(wp):: CL,bl1,bl2,r3l,time,Ek
      real(wp),allocatable:: massi(:,:)
      integer:: j,i,m,n,ib,k,iat,ii,imve,ic=0,natomc,ix,iy,iz,istep
      integer:: n_crossbond_x,n_crossbond_y,n_crossbond_z,n_silicon_etch,itemp
      logical:: connected_x,connected_y,connected_z,ok,spanning_cluster
      type(atom_list):: lst_positive_x,lst_negative_x,lstSi
      type(atom_list):: lst_positive_y,lst_negative_y
      type(atom_list):: lst_positive_z,lst_negative_z
	  real(wp),allocatable:: coprxyz(:,:)
      character(32):: ctmp
!
      open(unit=16,file='connect_struct.out')
      open(unit=14,file='simbox.xyz')
      open(unit=15,file='neighbor.out')
      open(unit=17,file='energy.out')
      n_silicon_etch = 320
      natom_max = natom_unit*nunitc*nunitc*nunitc
      nseed = 0
      allocate(crossbond_x(natom_max),crossbond_y(natom_max),crossbond_z(natom_max))
      allocate(rxyz(1:natom_max,3),coprxyz(natom_max,3))
      allocate(fxyz(1:natom_max,3))
      allocate(vxyz(1:natom_max,3))
      allocate(atom(1:natom_max))
      allocate(proximity(natom_max,4))
      allocate(massi(natom_max,3))
      allocate(rxyzc(1:natom_unit,3))
      proximity = 0
      write (*,*) natom_max
      CL = 0.7168_wp*0.162_wp/0.1552_wp

! read input data
      open(unit=13,FILE='unitcell.dat')
      do i=1,natom_unit
         read(13,*) ctmp,rxyz(i,1),rxyz(i,2),rxyz(i,3)
         atom(i) = name2atom(ctmp)
         if(atom(i) == iSilicon) then
            massi(i,:) = 1.0_wp/mSi 
         else
            massi(i,:) = 1.0_wp/mOx
         end if
      end do
     rxyz=rxyz/10.0_wp
      close(13,status='keep')

!------------------------------------------------------------
      boxl = nunitc*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      bl1 = 0.176_wp
      T(1) = 15000.0_wp
      T(2) = 10000.0_wp
      T(3) = 9000.0_wp
      T(4) = 8000.0_wp
      T(5) = 7000.0_wp
      T(6) = 5000.0_wp
      T(7) = 3000.0_wp
      T(8) = 1000.0_wp
      T(9) = 500.0_wp
      T(10) = 300.0_wp
!
! Creating the simulation box
!
      natom = 0
      rxyzc(1:natom_unit,:) = rxyz(1:natom_unit,:)
      do k = 1,nunitc
      do j = 1,nunitc
      do i = 1,nunitc
      do m=1,natom_unit
         natom = natom + 1
         atom(natom) = atom(m)
           
         if (atom(natom) == iSilicon) then
            massi(natom,:) = 1.0_wp/mSi
         else 
            massi(natom,:) = 1.0_wp/mOx
         end if        
          
         rxyz(natom,1) = rxyzc(m,1)+(i-1)*CL
         rxyz(natom,2) = rxyzc(m,2)+(j-1)*CL
         rxyz(natom,3) = rxyzc(m,3)+(k-1)*CL
      end do
      end do
      end do
      end do
      write (*,*) natom
! origin at centre of box
      rxyz(1:natom,1:3)=rxyz(1:natom,1:3) - CL*0.5_wp*(nunitc-1)
      forall(i=1:natom) rxyz(i,1:3)=rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)
!      call write_xyz(14,0)
!      close(14,status='KEEP')

! determine connectivity
      do k=1,natom
         n=0
         do j=1,natom
            if(j == k) cycle
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            r3(:)=r3(:)-boxl*anint(r3(:)*boxli)
            if ( (dot_product(r3,r3) <= bl1*bl1) ) then
                n=n+1
                proximity(k,n)=j
            end if
         end do
      end do
      do i=1,natom
         write (15,*) atom_name(atom(i)),proximity(i,1:4)
      end do


      call INIT_NLIST(boxl,boxl,boxl,3.0_wp/angstrom)
      call NEW_NLIST

      close(15,status='keep')   
      call COMVEL(natom,T(1),vxyz)
      vxyz = vxyz*sqrt(massi)
      call set_bond_angle_lists()
      call force_keating
      call repul_energy(1,natom)

      time = 0.0_wp
      istep = 0
      dtreal = 0.01_wp
      dt = dtreal*dt_fs
      do itemp =2,10
         vxyz(1:natom,:)=vxyz(1:natom,:)*sqrt(T(itemp)/T(itemp-1))
         do n=1, number_moves
            istep = istep + 1
               coprxyz = rxyz
            call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
            forall(i=1:natom) rxyz(i,1:3)=rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)

            ! check if the neighbour list needs to be updated
                  do i = 1,natom
               IF (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3)))then
                  call NEW_NLIST
                  EXIT
               END IF
                  end do

            call force_keating
            call repul_energy(1,natom)
          
            call vverlet_b(dt,massi,vxyz,fxyz)         
            time = time + dtreal
            if(mod(istep,100)==0) then
!  print *,sum(fxyz(1:natom,:),dim=1)
              Ek = 0.0_wp
              do m=1, natom
              Ek = Ek + (0.5_wp/massi(m,1))*dot_product(vxyz(m,:),vxyz(m,:))
              end do
            write (17,*) time, Ek+energy+repulenergy/2.0_wp, T(itemp)
            call flush(17)
            end if 
         end do
      end do
      close(17,status='keep')

      call write_xyz(14,0)
      close(14,status='KEEP')

end program simbox




