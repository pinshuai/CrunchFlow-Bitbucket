SUBROUTINE SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz)
USE crunchtype
USE concentration
USE medium
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: npot
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: npt

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: faraday
REAL(DP)                                                    :: gramsperL
REAL(DP)                                                    :: portemp
REAL(DP)                                                    :: convert
REAL(DP)                                                    :: correct

faraday = 96485.0       
surfcharge = 0.0

DO npt = 1,npot
  is = ispot(npt)
  k = ksurf(is)

  IF (volin(k,jinit(jx,jy,jz)) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN
    correct = wtmin(k)*specific(k,jinit(jx,jy,jz))*voltemp(k,jinit(jx,jy,jz))/volmol(k)   !!  m^2 mineral/m^3 BV
  ELSE
    correct = wtmin(k)*specific(k,jinit(jx,jy,jz))*volfx(k,jx,jy,jz)/volmol(k)   !!  m^2 mineral/m^3 BV
  END IF

  surfcharge(k) = surfcharge(k) + zsurf(is)*spsurf10(is,jx,jy,jz)*faraday/correct
  DO ns = 1,nsurf_sec
    IF (islink(ns) == is) THEN
      surfcharge(k) = surfcharge(k) + zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct    !!  Equation 2.1 in Dzombak
    END IF
  END DO

END DO

RETURN
END SUBROUTINE SurfaceCharge
