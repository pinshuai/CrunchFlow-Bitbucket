SUBROUTINE SurfaceCharge_init(ncomp,nspec,nsurf,nsurf_sec,npot,portemp,nco)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: npot
INTEGER(I4B), INTENT(IN)                                    :: nco

REAL(DP), INTENT(IN)                                        :: portemp

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: npt

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: faraday
REAL(DP)                                                    :: gramsperL
REAL(DP)                                                    :: term1


faraday = 96485.0       
surfcharge_init = 0.0

DO npt = 1,npot
  is = ispot(npt)
  k = ksurf(is) 

  IF (volin(k,nco) == 0.0d0) THEN
    gramsperL = wtmin(k)*voltemp(k,nco)/(rocond(nco)*volmol(k)*portemp)  ! units of g/L    
  ELSE
    gramsperL = wtmin(k)*volin(k,nco)/(rocond(nco)*volmol(k)*portemp)  ! units of g/L
  END IF

  term1 = faraday/(specific(k,nco)*gramsperL)

  surfcharge_init(k) = surfcharge_init(k) + zsurf(is)*spsurftmp10(is)*term1
  DO ns = 1,nsurf_sec
    IF (islink(ns) == is) THEN
      surfcharge_init(k) = surfcharge_init(k) + zsurf(ns+nsurf)*spsurftmp10(ns+nsurf)*term1
    END IF
  END DO
END DO

RETURN
END SUBROUTINE SurfaceCharge_init
