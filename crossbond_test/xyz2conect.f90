
include 'precision_mod.f90'
INCLUDE 'command_line_mod.f90'
INCLUDE 'command_arg_mod.f90'

MODULE global_vars_mod
   USE precision_mod, only: wp
!   USE phys_const_mod
   real(wp),parameter:: Angstrom = 1.0_wp  ! Angstroms/nm
   real(wp),parameter:: ulength = 1.0e-10_wp*Angstrom
!   real(wp),parameter:: uenergy = eV
!   real(wp),parameter:: upressure = uenergy/(ulength**3)
   real(wp),parameter:: torr2Pa = 101325.0_wp/760.0_wp
   real(wp),parameter:: pi = 3.1415926535897932384626433832795029_wp
   real(wp),parameter:: erg_ev = 6.241457E+11_wp
   real(wp),parameter:: K_ev = 8.6173423E-5_wp
   real(wp),parameter:: qstar = 1.19999_wp
   real(wp):: del_rxn, e_activ, kBoltzT
   integer:: natom,natom_max
   integer:: n_GAS=0,n_atom_GAS,n_O2=0,n_atom_o2,n_N2=0,n_atom_n2
   integer:: n_CO2=0,n_atom_co2,n_h2o=0,n_atom_h2o
END MODULE global_vars_mod

include 'sort_mod.f90'
include 'files_mod.f90'
include 'seaton_mod.f90'
include 'atom_types_mod.f90'
INCLUDE 'Keating_parameters_mod_vonAlfthan.f90'
include 'constants_mod.f90'
include 'coordinates_mod.f90'
include 'connectivity_mod.f90'
INCLUDE 'check_structure_mod_Silica.f90'
include 'nlist_mod3.f90'
include 'bond_angle_types_mod.f90'
include 'bond_angle_list_mod.f90'
include 'frames_mod.f90'
include 'readline_mod.f90'
include 'rand_mod.f90'
include 'matvec3d_mod.f90'
include 'HKNonLattice2.f90'
include 'list_mod.f90'


PROGRAM crossbond_correct
      USE precision_mod, only: wp
      USE command_line_mod
      USE command_arg_mod
      USE global_vars_mod
      USE bond_angle_list_mod
      USE constants_mod
      USE atom_types_mod
      USE connectivity_mod
      USE coordinates_mod
      USE files_mod
      USE seaton_mod
      USE readline_mod
      USE nlist_mod
      USE frames_mod
      USE HKNonLattice_mod
      USE list_mod
      USE check_structure_mod
      IMPLICIT NONE
      integer,parameter:: nl = 1 
      integer,parameter:: nunitc(3) = 2*nl + 1
      integer:: ip,itmp,ifirst
      integer:: i,j,jj,k,kk,m,n,iox,io_config,io_topol,narg,stat
      integer:: iu,ios,ia,ig,ilast,iatom,natom_orig,ix,iy,iz,ib
      integer:: n_crossbond,ncross,nc,ic,neigh2,neighmx=125
      integer,allocatable:: atomc(:),proximity_temp(:,:)
      integer,allocatable:: cross_atom(:),atom_temp(:),crossbond(:)
      real(wp):: timep(0:10),bond_lens(100)
      real(wp):: r3(3),rj(3),dr(3),bondl_orig=1.55_wp
      real(wp),allocatable:: rxyz_temp(:,:)
      real(wp):: CL,CL2,x,y,z,en
      character(len= 132):: coord_file,con_file
      character(len=32):: carg*132,word,ctmp
      integer:: silst(100),nsil,j1,j2,ia1(1),ia2(1),nn,neighbors(125)
      logical:: consistent
      narg = command_argument_count()
      call next_command_argument(carg, "command",stat)
      if (narg < 1) then
         write(*,*)'usage :'
         write(*,*) trim(carg),' atom_file'
         stop
      end if
      call next_command_argument(coord_file, "input file",stat)
      nc = len_trim(coord_file)
      con_file = coord_file(1:nc-3)//'con'  ;print *,'con_file = ',trim(con_file)
      call get_free_file_unit(io_config)
      open(unit=io_config,file=trim(coord_file))
      call get_free_file_unit(io_topol)
      open(unit=io_topol,file=trim(con_file))

!
      read (io_config,*) natom
      read (io_config,*) boxl
!
      boxl2 = boxl/2.0_wp
      boxli = 1.0_wp/boxl

      natom_max = natom
      allocate(rxyz(natom_max,3))
      allocate(atom(natom_max))
      atom = 0
      allocate(proximity(natom_max,4))
      proximity = 0
      do i = 1,natom
         read (io_config,*) ctmp,(rxyz(i,j),j = 1,3)
         atom(i) = name2atom(trim(ctmp))
      end do
!

      CALL cpu_time(timep(0))

      call INIT_NLIST(-boxl2,boxl2,5.0_wp,natom_max)
      call NEW_NLIST(1,natom)
!
! Loop over the Oxygen atoms
      oxygen_loop: do iox = 1,natom
         if (atom(iox) /= iOxygen) cycle oxygen_loop
         r3(1:3) = rxyz(iox,1:3)
         ic = CELL(r3)
         CALL NEIGCELL(ic,1,nn,neighbors)
         nsil = 0
         silst = 0
         bond_lens = 0.0_wp
         cell_loop: do jj = 1, nn
            nc = neighbors(jj)
            if(nc==0) cycle cell_loop
            j=HOC(nc)
            cell_atom_loop: do while (j/=0)   ! loop over cell Silicon atoms
               if(atom(j) /= iSilicon) GOTO 100
               if (j==iox) GOTO 100
if(ALL(proximity(j,:)/=0)) GOTO 100
               dr = r3-rxyz(j,1:3)
               CALL pbc(dr)
               if (dot_product(dr,dr) < (1.5*bondl_SiO)**2) then
                  nsil=nsil+1
                  bond_lens(nsil) = dot_product(dr,dr)
                  silst(nsil) = j 
               end if
100            CONTINUE
               j=LL(j)
            end do cell_atom_loop
         end do cell_loop

         if (nsil < 2) then
            print*, 'defect zone,nsil = ', nsil
            print*,'iox = ',iox
            !stop
         end if
         ia1 = minloc(bond_lens, (bond_lens > 0.0_wp))
if(ia1(1)>0)         ia2 = minloc(bond_lens, (bond_lens > bond_lens(ia1(1))))
if(ia1(1)>0)         j1 = silst(ia1(1))
if(ia2(1)>0)         j2 = silst(ia2(1))

!print *,'iox = ',iox
!do i = 1,nsil
!print *,i,sqrt(bond_lens(i)),silst(i)
!end do

         if(j1 > j2)then
            jj = j1
            j1 = j2
            j2 = jj
         end if

!print *,ia1,ia2
!print *,j1,j2
!print *,'prox iox = ',proximity(iox,:)
!print *,'prox j1  = ',proximity(j1,:)
!print *,'prox j2  = ',proximity(j2,:)

if(ia1(1)>0)         call set_proximity(iox,0,j1)
if(ia1(1)>0)         call set_proximity(j1,0,iox)
if(ia2(1)>0)         call set_proximity(iox,0,j2)
if(ia2(1)>0)         call set_proximity(j2,0,iox)
!         if (nsil > 2) then 
!            print*, 'defect zone,nsil = ', nsil
!         end if
      end do oxygen_loop
      
      CALL cpu_time(timep(1))
      print*, 'time taken = ', timep(1) - timep(0)
      
      call check_proximity(1,natom,consistent,ifirst)
      if (.NOT. consistent) then
         print*, 'not consistent'
         call write_atom_info(6,ifirst)
      end if
      call check_proximity2(1,natom,atom,proximity,consistent,ifirst)
      if (.NOT. consistent) then
         print*, 'not consistent 2'
         call write_atom_info(6,ifirst)
      end if

      do i = 1,natom
         select case(atom(i))
         case(iSilicon)
            if(count(proximity(i,:)>0) /= 4)then
               print *,'proximity error'
               print *,'i = ',i
               print *,(proximity(i,j),j=1,4)
            end if
         case(iOxygen)
            if(count(proximity(i,:)>0) /= 2)then
               print *,'proximity error'
               print *,'i = ',i
               print *,(proximity(i,j),j=1,4)
            end if
            if(is_2si_ring(i))then
               print *,i,' is in a 2-si ring'
               print *,i,'   ',(proximity(i,j),j=1,2)
               print *,proximity(i,1),'   ',(proximity(proximity(i,1),j),j=1,4)
               print *,proximity(i,2),'   ',(proximity(proximity(i,2),j),j=1,4)
            end if
         end select
      end do

      call init_hknonlattice(natom)
      call hknonlattice(natom,proximity,n_cluster,atoml)
      print *,'n_cluster = ',n_cluster
      if(n_cluster /= 1)then
         stop 'there should be one cluster'
      end if
      
      write(io_topol,*) natom
      write(io_topol,*) 'conectivity'
      do i = 1,natom
         write(io_topol,"(5(1x,i6))") i,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
END PROGRAM crossbond_correct

