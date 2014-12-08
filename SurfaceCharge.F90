SUBROUTINE SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz,time)
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
REAL(DP), INTENT(IN)                                        :: time

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
REAL(DP)                                                    :: volMinimum

faraday = 96485.0d0     
surfcharge = 0.0d0

DO npt = 1,npot

  k = kpot(npt)

  IF (volin(k,jinit(jx,jy,jz)) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN
    correct = wtmin(k)*specific(k,jinit(jx,jy,jz))*voltemp(k,jinit(jx,jy,jz))/volmol(k)   !!  m^2 mineral/m^3 BV
  ELSE
    volMinimum = volfx(k,jx,jy,jz)
    if (volMinimum < 1.0D-15) then
       volMinimum = 1.0D-15
    end if
    correct = wtmin(k)*specific(k,jinit(jx,jy,jz))*volMinimum/volmol(k)   !!  m^2 mineral/m^3 BV
  END IF

  if (correct == 0.0) then
    write(*,*) ' Divide by zero in surface charge'
    read(*,*)
  end if

  DO is = 1,nsurf
    IF ( ksurf(is) == kpot(npt) ) THEN
      surfcharge(k) = surfcharge(k) + zsurf(is)*spsurf10(is,jx,jy,jz)*faraday/correct
    END IF
  END DO

  DO ns = 1,nsurf_sec
    IF (ksurf(islink(ns)) == kpot(npt)) THEN
      surfcharge(k) = surfcharge(k) + zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct    !!  Equation 2.1 in Dzombak
    END IF
  END DO


END DO


RETURN
END SUBROUTINE SurfaceCharge
