SUBROUTINE jacsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE
!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz

!  Internal variables

REAL(DP)                                             :: sum
REAL(DP)                                             :: mutemp
REAL(DP)                                             :: surftemp

INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: iss
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: is

fsurf_local = 0.0

DO ns = 1,nsurf_sec
  surftemp = spsurf10(ns+nsurf,jx,jy,jz)
  DO is = 1,nsurf
    iss = ncomp + is
    IF (musurf(ns,iss) /= 0.0) THEN
      mutemp = musurf(ns,iss)
      DO i2 = 1,ncomp+nsurf
        fsurf_local(is,i2) = fsurf_local(is,i2) +    &
           mutemp*musurf(ns,i2)*surftemp
      END DO
    END IF
  END DO
END DO

DO is = 1,nsurf
  iss = ncomp + is
  fsurf_local(is,iss) = fsurf_local(is,iss) + spsurf10(is,jx,jy,jz)
END DO

RETURN
END SUBROUTINE jacsurf_local
!***********************************************************
