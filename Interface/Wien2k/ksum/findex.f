function FIndex(jatom, iatom, natom)
  IMPLICIT NONE
  INTEGER :: FIndex
  INTEGER, intent(in) :: iatom(natom)
  INTEGER, intent(in) :: jatom, natom
  ! locals
  INTEGER :: i
  FIndex=0
  do i=1,natom
     if (iatom(i).eq.jatom) then
        FIndex=i
        return
     endif
  enddo
  return
end function FIndex
