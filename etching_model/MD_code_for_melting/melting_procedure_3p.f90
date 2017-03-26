
include 'precision_mod.f90'
include 'files_mod.f90'
include 'rand_mod.f90'
include 'sort_mod.f90'
include 'global_vars.f90'
include 'HKNonLattice2.f90'
include 'seaton_mod.f90'
include 'Keating_parameters_vonAlfthan.f90'
include 'constants_mod.f90'
include 'atom_types.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
include 'atom_list_mod.f90'
include 'bond_list_mod.f90'
include 'frames_mod.f90'
include 'comvel_mod.f90'
include 'force_keating.f90'
!include 'force_angle_mm.f90'
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
      logical,parameter:: make_movie = .TRUE.
      integer,parameter:: natom_unit = 3
      integer,parameter:: nunitc = 1
      integer,parameter:: number_moves = 2000000
      real(wp),parameter:: mSi = 28.06_wp
      real(wp),parameter:: mOx = 16.00_wp
      real(wp):: r3(3),T(0:10),dt,dtreal
      real,allocatable:: rxyzc(:,:)
      real(wp):: CL,bl1,Ek
      real(wp),allocatable:: massi(:,:)
      integer:: j,i,ia,n,k,it,imve,istep,itemp,ntemp,narg
      real(wp),allocatable:: coprxyz(:,:)
      character(32):: ctmp
!
      narg = command_argument_count()
      if(narg < 1)then
         write(*,*)'usage :'
         write(*,*)'  melting_procedure timestep_fs'
         stop
      end if
      call get_command_argument (1, ctmp)
      read(unit=ctmp,fmt=*) dtreal
      open(unit=16,file='connect_struct.out')
      open(unit=14,file='simbox.xyz')
      open(unit=15,file='neighbor.out')
      open(unit=17,file='energy.out')
!
      natom_max = natom_unit*nunitc*nunitc*nunitc
      nseed = 0
      time_s = length_m*sqrt(amu/ev)
      dt_fs = 1.0e-15/time_s
      dtreal = 0.01_wp ! in femtoseconds
      dt = dtreal*dt_fs
!
print *,'unit of time = ',time_s
print *,'dt_fs = ',dt_fs
print *,'dt = ',dtreal,' fs'
print *,'dt = ',dt
print *,'dt*unit of time = ',dt*time_s

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
!      do i=1,natom_unit
!         read(13,*) ctmp,rxyz(i,1),rxyz(i,2),rxyz(i,3)
!         atom(i) = name2atom(ctmp)
!         if(atom(i) == iSilicon) then
!            massi(i,:) = 1.0_wp/mSi 
!         else
!            massi(i,:) = 1.0_wp/mOx
!         end if
!      end do
!      rxyz(1:natom_unit,:) = rxyz(1:natom_unit,:)/10.0_wp
      close(13,status='keep')
      rxyz(1,:) = (/0.0_wp,bondl_SiO,0.0_wp /)
      rxyz(2,:) = 0.0_wp
      rxyz(3,:) = (/bondl_SiO,0.0_wp,0.0_wp /)
      atom(1:3) = (/iOxygen,iSilicon,iOxygen /)
      massi(1,:) = 1.0_wp/mOx
      massi(2,:) = 1.0_wp/mSi
      massi(3,:) = 1.0_wp/mOx
      vxyz(1,:) = (/  0.01_wp, 0.01_wp,0.0_wp /)
      vxyz(3,:) = (/  0.01_wp, 0.01_wp,0.0_wp /)
      vxyz(2,:) = -2.0_wp*(mOx/mSi)*vxyz(1,:)
      natom = 3
      vxyz = 0.0_wp
      
!------------------------------------------------------------
      boxl = nunitc*CL
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl
      bl1 = 0.176_wp
      ntemp = 1
      T(0) = 12.0_wp
      T(1) = 10.0_wp
      T(2) =  9000.0_wp
      T(3) =  8000.0_wp
      T(4) =  7000.0_wp
      T(5) =  6000.0_wp
      T(6) =  5000.0_wp
      T(7) =  3000.0_wp
      T(8) =  1000.0_wp
      T(9) =   500.0_wp
      T(10) =  300.0_wp
!
! Creating the simulation box
!
!      natom = 0
!      rxyzc(1:natom_unit,:) = rxyz(1:natom_unit,:)
!      do k = 1,nunitc
!      do j = 1,nunitc
!      do i = 1,nunitc
!      do ia = 1,natom_unit
!         natom = natom + 1
!         atom(natom) = atom(ia)
!         if (atom(natom) == iSilicon) then
!            massi(natom,:) = 1.0_wp/mSi
!         else 
!            massi(natom,:) = 1.0_wp/mOx
!         end if         
!         rxyz(natom,1) = rxyzc(ia,1)+(i-1)*CL
!         rxyz(natom,2) = rxyzc(ia,2)+(j-1)*CL
!         rxyz(natom,3) = rxyzc(ia,3)+(k-1)*CL
!      end do
!      end do
!      end do
!      end do
      write (*,*) natom
! origin at centre of box
      rxyz(1:natom,1:3)=rxyz(1:natom,1:3) - CL*0.5_wp*(nunitc-1)
      forall(i=1:natom) rxyz(i,1:3)=rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)
!     call write_xyz(14,0)
!     close(14,status='KEEP')

! determine connectivity
      do k = 1,natom
         n = 0
         do j=1,natom
            if(j == k) cycle
            if(atom(k) == atom(j)) cycle
            r3(:)=rxyz(j,:)-rxyz(k,:)
            r3(:)=r3(:)-boxl*anint(r3(:)*boxli)
            if ( (dot_product(r3,r3) <= bl1*bl1) ) then
                n = n+1
                proximity(k,n)=j
            end if
         end do
      end do
      do i=1,natom
         write (*,*) atom_name(atom(i)),proximity(i,1:4)
      end do
!     close(15,status='keep') 

! initilize the neighbour list
      call INIT_NLIST(boxl,boxl,boxl,3.0_wp/angstrom)
      call NEW_NLIST

  
      call COMVEL(natom,T(1),vxyz)
      vxyz = vxyz*sqrt(massi)
      call set_bond_angle_lists

      if (make_movie) then ! write out first frame of movie
         call new_frame_file(imve,'frame',0)
         call write_frame(imve,1-nseed,natom)
      end if

      
      fxyz(1:natom,:) = 0.0_wp
      call force_keating_bond
      call force_keating_angle
!      call repul_energy(1,natom)

      istep = 0
      temperature: do itemp =1,ntemp

         vxyz(1:natom,:) = vxyz(1:natom,:)*sqrt(T(itemp)/T(itemp-1))

         md_loop: do it = 1, number_moves
            istep = istep + 1
            coprxyz(1:natom,:) = rxyz(1:natom,:)

            call vverlet_a(dt,massi,rxyz,vxyz,fxyz)
            forall(i=1:natom) rxyz(i,1:3) = rxyz(i,1:3)-boxl*anint(rxyz(i,1:3)*boxli)

            ! check if the neighbour list needs to be updated
            do i = 1,natom
               if (CELL(rxyz(i,1:3)) /= CELL(coprxyz(i,1:3)))then
                  call NEW_NLIST
                  exit
               end if
            end do

            fxyz(1:natom,:) = 0.0_wp
            call force_keating_bond
            call force_keating_angle
!            call repul_energy(1,natom)
          
            call vverlet_b(dt,massi,vxyz,fxyz)         

            if(mod(istep,50)==0) then
               print *,sum(fxyz(1:natom,:),dim=1)
               Ek = 0.0_wp
               do i = 1,natom
                  Ek = Ek + (0.5_wp/massi(i,1))*dot_product(vxyz(i,:),vxyz(i,:))
               end do
               write(17,'(6g18.10)') istep*dtreal, Ek+ebond+eangle+repulenergy*0.5_wp, &
                            ek,ebond+eangle,repulenergy*0.5_wp,T(itemp)
               call flush(17)
            end if
            if (make_movie.and.mod(istep,200)==0.and. (istep/200 < 9999) ) then
               call new_frame_file(imve,'frame',istep/200)
               call write_frame(imve,1,natom)
            end if

         end do md_loop
         write (17,'(/)')

      end do temperature

      close(17,status='keep')
      call write_xyz(14,0)
      close(14,status='KEEP')

end program simbox




