
MODULE select_atom_mod
   USE precision_mod
   USE global_vars
   USE constants_mod
   USE atom_types
   USE coordinates_mod
   USE connectivity_mod
   USE rand_mod
   implicit none
CONTAINS

   SUBROUTINE get_atom_for_attachment4(iatom)
      integer,intent(out):: iatom
      integer:: arraybonds(natom_max),n,L
      n = 0
      arraybonds(1:natom) = 0
      do L = 1,natom
         if (atom(L) /= iOxygenH) cycle
         n = n + 1
         arraybonds(n) = L
      end do
      iatom = arraybonds( int(rand()*n)+1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment3(iz,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      integer,intent(out):: iatom
      integer:: arraybonds(natom_max),n,L,ic,nc
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            n = n + 1
            arraybonds(n) = L
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      iatom = arraybonds( int(rand()*n)+1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment2(iz,diam,bondl,iatom)
      USE nlist_mod
      integer,intent(in):: iz
      real(wp),intent(in):: diam,bondl
      integer,intent(out):: iatom
      integer,parameter:: ntrial_max = 100
      real(wp):: r3v(3),alfa,beta
      integer:: arraybonds(natom_max),i,n,L,it,ic,nc
      n = 0
      arraybonds(1:natom) = 0
      call Z_NEIGH_CELL(iz,neigh,ncell)
      cell_loop: do ic = 1,neigh
         nc = ncell(ic)
         if (nc == 0)cycle cell_loop
         L = HOC(nc)
         atom_loop: do while (L /= 0)
            if (atom(L) /= iOxygenH) GOTO 100
            trial_loop: do it = 1,ntrial_max
               alfa = acos(2.0_wp*rand() - 1.0_wp)
               beta = rand()*2.0_wp*pi
               rxyz(natom+1,3) = rxyz(L,3) + bondl*cos(beta)
               rxyz(natom+1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
               rxyz(natom+1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
               call pbc(rxyz(natom+1,:))
               do i = 1,natom
                  if (atom(i) == iHydrogen) cycle
                  if (atom(i) == iSilicon) cycle
                  if (i == L) cycle
                  r3v(1:3) = rxyz(i,1:3) - rxyz(natom+1,1:3)
                  call pbc(r3v)
                  if (dot_product(r3v,r3v) <= diam**2) CYCLE trial_loop
               end do
               exit trial_loop
            end do trial_loop
            if (it < ntrial_max) then
               n = n + 1
               arraybonds(n) = L
            end if
100         L = LL(L)
         end do atom_loop
      end do cell_loop
      iatom = arraybonds( int(rand()*n)+1 )
   END SUBROUTINE

   SUBROUTINE get_atom_for_attachment(z_dep,diam,bondl,iatom)
      real(wp),intent(in):: z_dep,diam,bondl
      integer,intent(out):: iatom
      integer,parameter:: ntrial_max = 100
      real(wp):: r3v(3),alfa,beta
      integer:: arraybonds(natom_max),i,n,L,it,iv(1)
!
      n = 0
      arraybonds(1:natom) = 0
      atom_loop: do L = 1,natom
         if (atom(L) /= iOxygenH) cycle atom_loop
         trial_loop: do it = 1,ntrial_max
            alfa = acos(2.0_wp*rand() - 1.0_wp)
            beta = rand()*2.0_wp*pi
            rxyz(natom+1,3) = rxyz(L,3) + bondl*cos(beta)
            rxyz(natom+1,1) = rxyz(L,1) + bondl*cos(alfa)*sin(beta)
            rxyz(natom+1,2) = rxyz(L,2) + bondl*sin(alfa)*sin(beta)
            call pbc(rxyz(natom+1,:))
            do i = 1,natom
               if (atom(i) == iHydrogen) cycle
               if (atom(i) == iSilicon) cycle
               if (i == L) cycle
               r3v(1:3) = rxyz(i,1:3) - rxyz(natom+1,1:3)
               call pbc(r3v)
               if (dot_product(r3v,r3v) <= diam**2) CYCLE trial_loop
            end do
            exit trial_loop
         end do trial_loop
         if (it >= ntrial_max)cycle atom_loop
         n = n + 1
         arraybonds(n) = L
      end do atom_loop
      iv = MINLOC(ABS(rxyz(arraybonds(1:n),3)-z_dep))
      iatom = arraybonds(iv(1))
      if (iatom <= nseed) then
         do
            iatom = rand()*nseed + 1.0
            if (atom(iatom) == iOxygenH) exit
         end do
      end if
   END SUBROUTINE get_atom_for_attachment

END MODULE select_atom_mod

