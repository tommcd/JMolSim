
MODULE atom_types
   USE precision_mod, only: wp
   implicit none
   integer,parameter:: iOxygen = 0
   integer,parameter:: iSilicon = 1
   integer,parameter:: iOxygenH = 2
   integer,parameter:: iHydrogen = 3
   integer,parameter:: iOw=4
   integer,parameter:: iHw=5
   integer,parameter:: ntyp = 5
   character(2),parameter:: atom_name(0:ntyp)=(/ 'O ','Si','OH','H ','Ow','Hw' /)
   integer,parameter:: ncmax(0:ntyp)=(/2,4,2,1,1,2/)
   integer,allocatable:: atom(:),Ox_atom(:,:),N2_atom(:,:),CO2_atom(:,:)
CONTAINS

   pure function name2atom(c)
      integer:: name2atom
      character(*),intent(in):: c
      integer:: i
      do i = 0,ntyp
         if(trim(c) == trim(atom_name(i)))then
            name2atom = i
            exit
         end if
      end do
      if (i > ntyp) name2atom = -1
   end function

END MODULE

