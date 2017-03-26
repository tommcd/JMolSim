MODULE FRAMES_MOD
   USE precision_mod, only: wp
   implicit none
   public:: new_frame_file, write_frame

CONTAINS

   SUBROUTINE new_frame_file(iu,base_name,iframe)
      USE files_mod
      integer,intent(in):: iframe
      integer,intent(out):: iu
      character(len=*),intent(in):: base_name
      character(len=32):: ctmp
      character(len=132):: frame_file
      integer:: ios
      logical:: lexists,lopen
      call get_free_file_unit(iu)
      write(unit=ctmp,fmt='(i4.4)') iframe
      frame_file = trim(base_name)//'_'//trim(adjustl(ctmp))//'.pdb'
      inquire(file=trim(frame_file),iostat=ios,exist=lexists,opened=lopen)
      if (ios == 0) then
         if (lexists .and. lopen) close(iu)
      else
         write(*,*)' new_frame_file: error, iostat = ',ios
         stop
      end if
      open(unit=iu,file=trim(adjustl(frame_file)),form='formatted')
   END SUBROUTINE new_frame_file

   SUBROUTINE write_frame(iu,ifirst,ilast)
      USE files_mod
      USE bond_angle_list_mod
      USE global_vars_mod, only: angstrom
      USE connectivity_mod
      USE atom_types_mod, only: atom,atom_name,ncmax
      USE coordinates_mod, only: rxyz,boxl2,boxl
      integer,intent(in):: iu,ifirst,ilast
      integer:: i,j,ib,ia,iu2
      logical:: lopen,fix_bonds
      character(len=132):: frame_file
      write(iu,'(a6,2x,i7)') 'REMARK',ilast - ifirst + 1
      write(iu,'(a6,1x,3(1x,f8.4))') 'REMARK',boxl*Angstrom
      do i = ifirst,ilast
         write(iu,110)'HETATM',i,atom_name(atom(i)),'   ',i,rxyz(i,1:3)*angstrom
      end do
      call bond_list
!
      inquire(unit=iu,name=frame_file)
      j = len_trim(frame_file)
      frame_file = frame_file(1:j-3)//'scr'
      fix_bonds = .false.
      do i = 1,nbondtot
         ia = ibond(1,i)
         ib = ibond(2,i)
         if (ANY( ABS(rxyz(ia,1:3) - rxyz(ib,1:3)) > 10.0_wp/angstrom )) then
            fix_bonds = .true.
            inquire(file=frame_file,opened=lopen)
            if (.not.lopen) then
               call get_free_file_unit(iu2)
               open(unit=iu2,file=trim(frame_file))
            end if
            write(iu2,*)'select ',ia,',',ib
            write(iu2,*)'wireframe off'
         end if
!         write(iu,"(a6,5i5)")'CONECT',ibond(1,i),ibond(2,i)
      end do
      if (fix_bonds) then
         write(iu2,*)
         close(iu2,status='keep')
      end if
      do i = 1,ilast
         if ( all(proximity(i,:) == 0) ) cycle
         write(iu,"(a6,5i5)")'CONECT',i,(proximity(i,j),j = 1,ncmax(atom(i)))
      end do
      write(iu,110)'HETATM',99991,'AR','   ',99991,-boxl2(1)*angstrom,-boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99992,'AR','   ',99992,-boxl2(1)*angstrom, boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99993,'AR','   ',99993, boxl2(1)*angstrom,-boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99994,'AR','   ',99994, boxl2(1)*angstrom, boxl2(2)*angstrom,0.0
      write(iu,110)'HETATM',99995,'AR','   ',99995,-boxl2(1)*angstrom,-boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99996,'AR','   ',99996,-boxl2(1)*angstrom, boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99997,'AR','   ',99997, boxl2(1)*angstrom,-boxl2(2)*angstrom,5.0
      write(iu,110)'HETATM',99998,'AR','   ',99998, boxl2(1)*angstrom, boxl2(2)*angstrom,5.0
      write(iu,"(a6,5i5)")'CONECT',99991,99992
      write(iu,"(a6,5i5)")'CONECT',99991,99993
      write(iu,"(a6,5i5)")'CONECT',99991,99995
      write(iu,"(a6,5i5)")'CONECT',99992,99994
      write(iu,"(a6,5i5)")'CONECT',99992,99996
      write(iu,"(a6,5i5)")'CONECT',99993,99994
      write(iu,"(a6,5i5)")'CONECT',99993,99997
      write(iu,"(a6,5i5)")'CONECT',99994,99998
      write(iu,"(a6,5i5)")'CONECT',99995,99996
      write(iu,"(a6,5i5)")'CONECT',99995,99997
      write(iu,"(a6,5i5)")'CONECT',99996,99998
      write(iu,"(a6,5i5)")'CONECT',99997,99998
      write(iu,'(a/)')'END'
      close(iu,status='keep')
110   format(a6,i5,1x,a4,1x,a3,i6,4x,3f8.3)
   END SUBROUTINE write_frame

END MODULE FRAMES_MOD

