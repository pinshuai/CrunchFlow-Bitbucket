SUBROUTINE mindecay(ndecay,jx,jy,jz,delt)
USE crunchtype
USE params
USE concentration
USE mineral

IMPLICIT NONE

!  External arrays and variables

REAL(DP), INTENT(IN)                                       :: delt

INTEGER(I4B), INTENT(IN)                                   :: ndecay
INTEGER(I4B), INTENT(IN)                                   :: jx
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz

!  Internal variables and arrays

REAL(DP)                                                   :: rdecay

INTEGER(I4B)                                               :: id
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: kd
INTEGER(I4B)                                               :: k
INTEGER(I4B)                                               :: isotope

!  NOTE:  In initialization (mintope here is the stoichiometric coefficient for isotope)


DO id = 1,ndecay
!  Point to the primary species number
  i = idecay(id)
!  Cycle through minerals
  DO kd = 1,nmindecay(id)
!  Point to the mineral number
    k = kdecay(kd,id)
    
!  Cycle through isotopes of a particular decay species
    
    DO isotope = 1,nisotope(id)
      rdecay = -t_half/half_life(isotope,id)
      ratio_isotope(isotope,kd,id,jx,jy,jz) = ratio_isotope(isotope,kd,id,jx,jy,jz)*  &
                   EXP(-rdecay*delt)
    END DO
    
  END DO
END DO


!  Note:  gives directly the change in the reaction:
!     PuO2(s) = mumin(:,:)*Pu+++   Used to be 1, now some fraction of 1
!     Moles of radioactive species in mineral
!       molerad = mumin(i,k)*volfx(k,j)*volmol(k)
!     Moles of specific isotope in mineral
!       moleisotope = mintope(isotope,i,k)*volfx(k,j)*volmol(k)
!       call get_isotope(returns moles of isotope)

END SUBROUTINE mindecay
