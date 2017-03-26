MODULE bond_angle_types_mod
    USE precision_mod
    USE atom_types_mod
    USE constants_mod
    implicit none

CONTAINS

   SUBROUTINE bond_type(ia1,ia2,kbond,abond)
      integer,intent(in):: ia1,ia2
      real(wp),intent(out):: kbond,abond
      if (compare_bond(atom(ia1),atom(ia2),iSilicon,iOxygen)) then
         kbond = kSiO
         abond = aSiO
      else if (compare_bond(atom(ia1),atom(ia2),iSilicon,iOxygenH)) then
         kbond = kSiO
         abond = aSiO
      else if (compare_bond(atom(ia1),atom(ia2),iSilicon,iSilicon)) then
         kbond = KSiSi
         abond = ASiSi
         stop 'Si - Si bond not allowed'
      else if (compare_bond(atom(ia1),atom(ia2),iOxygenH,iHydrogen)) then
         kbond = KOH
         abond = AOH
      else
         stop 'unknown bond'
      end if
   END SUBROUTINE bond_type


   SUBROUTINE angle_type(ia1,ia2,ia3,kang,ctheta)
      integer,intent(in):: ia1,ia2,ia3
      real(wp),intent(out):: kang,ctheta
      if ( compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygen, iSilicon,iOxygen) .or.&
           compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygenH,iSilicon,iOxygenH).or.&
           compare_angle(atom(ia1),atom(ia2),atom(ia3),iOxygen, iSilicon,iOxygenH) ) then
         kang = KOSiO
         ctheta = ctheta_OSiO
      else if (compare_angle(atom(ia1),atom(ia2),atom(ia3),iSilicon,iOxygenH,iHydrogen)) then
         kang = KSiOH
         ctheta = cos(angle_SiOH)
      else if (compare_angle(atom(ia1),atom(ia2),atom(ia3),iSilicon,iOxygen,iSilicon)) then
         kang = KSiOSi
         ctheta = ctheta_SiOSi
      else
         stop 'unknown angle'
      end if
   END SUBROUTINE angle_type


   logical pure FUNCTION compare_bond(i1,i2,typ1,typ2)
      integer,intent(in):: i1,i2,typ1,typ2
      compare_bond = (i1 == typ1.and.i2 == typ2).or. &
                     (i1 == typ2.and.i2 == typ1)
   END FUNCTION compare_bond


   logical pure FUNCTION compare_angle(i1,i2,i3,typ1,typ2,typ3)
      integer,intent(in):: i1,i2,i3,typ1,typ2,typ3
      compare_angle = (i1 == typ1.and.i2 == typ2.and.i3 == typ3).or. &
                      (i1 == typ3.and.i2 == typ2.and.i3 == typ1)
   END FUNCTION compare_angle

END MODULE bond_angle_types_mod

