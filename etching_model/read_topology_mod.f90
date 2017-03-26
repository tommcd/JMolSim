
MODULE read_topology_mod
   implicit none
   private
   public :: read_topology
!
CONTAINS

   SUBROUTINE read_topology(inu,iall_angles,nmob,nap,nmw,nms,iout)
      USE precision_mod, only: wp
      USE dreiding_mod, only: nty
      USE angle_mob_mod, only: nmang,imang,rmtheta,rmkang
      USE angle_prot_mod, only: npang,ipang,rptheta,rpkang
      USE bond_mob_mod, only: nmbond,imbond,imblen,rmblen,rmdsqk
      USE bond_prot_mod, only: npbond,ipbond,ipblen,rpblen,rpdsqk
      USE connect_mob_mod ,only: imconn,imconnx
      USE connect_prot_mod ,only: ipconn,ipconnx
      USE atom_types_mod, only: itypem,itypep
      USE torsion_prot_mod, only: nptor,iptor,ityptor,rpvtor,rpphi0,rpcphi0
      implicit none
      integer,intent(in):: inu,iall_angles,iout,nmob,nap,nmw,nms
      integer:: i,j,ii,ia1,ia2,ia3,ia4,ib,ic,id,itst
      integer:: nindepang(4),npanga2(nap),nmanga2((nmw*1 + nms*5)*3)
!
! read in connectivity and types for atoms in 1 protein
      read(inu,*) itst
      if (itst /= nap) then
        write(*,*)'error reading in connectivity table'
        write(*,"(a,i4,a,i4)")'nap = ',nap,' /= ',itst
        stop
      end if
      do i = 1,nap
        read(inu,*) ii,itypep(i),ipconnx(0,i),  &
              (ipconnx(j,i),j = 1,ipconnx(0,i))
      end do

! read in connectivity and types for atoms mobile molecules
      read(inu,*) itst
      if (itst /= nmob) then
        write(*,*)'error reading in connectivity table'
        write(*,"(a,i4,a,i4)")'nmob = ',nmob,' /= ',itst
        stop
      end if
      do i = 1,nmob
        read(inu,*) ii,itypem(i),imconnx(0,i),  &
              (imconnx(j,i),j = 1,imconnx(0,i))
      end do

      nindepang(1) = 0
      nindepang(2) = 1   ! an atom connected to 2 other atoms has 1 indep
      nindepang(3) = 3   ! angle associated with it   ... etc.
      nindepang(4) = 5   ! in general   nindepang(i) = max(0,2*i - 3)
! assess connectivity for polymer matrix
      npbond = 0
      npang = 0
      nptor = 0
      ipconn(1:nap,1:nap) = -1      ! initialise connectivity matrix
      npanga2(1:nap) = 0      ! init array of no.of of angles around atom a2
      do ia1 = 1,nap
        ipconn(ia1,ia1) = 1
          do ib = 1,ipconnx(0,ia1)
          ia2 = ipconnx(ib,ia1)
          if (ia1 < ia2) then
            npbond = npbond + 1
            ipbond(1,npbond) = ia1
            ipbond(2,npbond) = ia2
            ipconn(ia1,ia2) = 2
            ipconn(ia2,ia1) = 2
          end if
          do ic = 1,ipconnx(0,ia2)
            ia3 = ipconnx(ic,ia2)
            if (ia1 == ia3) cycle    ! if atom_1 == atom_3 skip ang & tor
            if (ia1 < ia3) then
              ipconn(ia1,ia3) = 3
              ipconn(ia3,ia1) = 3
              if (iall_angles /= 0) then
                npanga2(ia2) = npanga2(ia2) + 1
                npang = npang + 1
                ipang(1,npang) = ia1
                ipang(2,npang) = ia2
                ipang(3,npang) = ia3
              else if (npanga2(ia2) < nindepang(ipconnx(0,ia2))) then
                npanga2(ia2) = npanga2(ia2) + 1
                npang = npang + 1
                ipang(1,npang) = ia1
                ipang(2,npang) = ia2
                ipang(3,npang) = ia3
              end if
            end if
            do id = 1,ipconnx(0,ia3)
              ia4 = ipconnx(id,ia3)
              if ((ia1 < ia4).and.(ia2 /= ia4)) then
              nptor = nptor + 1
              iptor(1,nptor) = ia1
              iptor(2,nptor) = ia2
              iptor(3,nptor) = ia3
              iptor(4,nptor) = ia4
              ipconn(ia1,ia4) = 4
              ipconn(ia4,ia1) = 4
              end if
            end do
          end do
        end do
      end do

! Connectivity for Salt/Water system (no torsion angles involved)
      nmbond = 0
      nmang = 0
      imconn(:,:) = -1
      nmanga2(:) = 0
      do ia1 = 1,nmob
        imconn(ia1,ia1) = 1
        do ib = 1,imconnx(0,ia1)
          ia2 = imconnx(ib,ia1)
          if (ia1 < ia2) then
            nmbond = nmbond + 1
            imbond(1,nmbond) = ia1
            imbond(2,nmbond) = ia2
            imconn(ia1,ia2) = 2
            imconn(ia2,ia1) = 2
          end if
          do ic = 1,imconnx(0,ia2)
            ia3 = imconnx(ic,ia2)
            if (ia1 == ia3) cycle          ! if atom_1 == atom_3 skip ang
            if (nmanga2(ia2) < nindepang(imconnx(0,ia2))) then
              if (ia1 < ia3) then
                nmanga2(ia2) = nmanga2(ia2) + 1
                nmang = nmang + 1
                imang(1,nmang) = ia1
                imang(2,nmang) = ia2
                imang(3,nmang) = ia3
                imconn(ia1,ia3) = 3
              end if
            end if
          end do
        end do
      end do
!
      if (iout > 0) then
        write(iout,"(/a)")'nmanga2:'
        do i = 1,nmob
          write(iout,"(2i5)") i,nmanga2(i)
        end do
        write(iout,"(/a)")'npanga2:'
        do i = 1,nap
          write(iout,"(2i5)") i,npanga2(i)
        end do
        write(iout,*)
      end if
      RETURN
   END SUBROUTINE read_topology
END MODULE read_topology_mod

