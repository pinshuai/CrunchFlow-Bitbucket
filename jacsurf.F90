SUBROUTINE jacsurf(ncomp,nsurf,nsurf_sec,nx,ny,nz)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
INTEGER(I4B), INTENT(IN)                             :: nx
INTEGER(I4B), INTENT(IN)                             :: ny
INTEGER(I4B), INTENT(IN)                             :: nz

!  Internal variables

REAL(DP)                                             :: sum
REAL(DP)                                             :: mutemp
REAL(DP)                                             :: surftemp

INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: is
INTEGER(I4B)                                         :: jx
INTEGER(I4B)                                         :: jy
INTEGER(I4B)                                         :: jz

fsurf = 0.0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      DO ns = 1,nsurf_sec
        surftemp = spsurf10(ns+nsurf,jx,jy,jz)
        DO is = 1,nsurf
          IF (musurf(ns,is+ncomp) /= 0.0) THEN
            mutemp = musurf(ns,is+ncomp)
            DO i2 = 1,ncomp+nsurf
              fsurf(is,i2,jx,jy,jz) = fsurf(is,i2,jx,jy,jz) +    &
                  mutemp*musurf(ns,i2)*surftemp
            END DO
          END IF
        END DO
      END DO

      DO is = 1,nsurf
        fsurf(is,is+ncomp,jx,jy,jz) = fsurf(is,is+ncomp,jx,jy,jz) + spsurf10(is,jx,jy,jz)
      END DO

    END DO
  END DO
END DO

RETURN
END SUBROUTINE jacsurf
!***********************************************************
