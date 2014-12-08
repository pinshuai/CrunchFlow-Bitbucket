SUBROUTINE jacpotential(ncomp,nsurf,nsurf_sec,npot,nx,ny,nz)
USE crunchtype
USE concentration
USE solver
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                                :: ncomp
INTEGER(I4B), INTENT(IN)                                                :: nsurf
INTEGER(I4B), INTENT(IN)                                                :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                                :: npot
INTEGER(I4B), INTENT(IN)                                                :: nx
INTEGER(I4B), INTENT(IN)                                                :: ny
INTEGER(I4B), INTENT(IN)                                                :: nz

!  Internal variables and arrays

REAL(DP)                                                                :: delta_z
REAL(DP)                                                                :: surfconc
REAL(DP)                                                                :: mutemp
REAL(DP)                                                                :: sum

INTEGER(I4B)                                                            :: npt2
INTEGER(I4B)                                                            :: is
INTEGER(I4B)                                                            :: is2
INTEGER(I4B)                                                            :: ns
INTEGER(I4B)                                                            :: i
INTEGER(I4B)                                                            :: jx
INTEGER(I4B)                                                            :: jy
INTEGER(I4B)                                                            :: jz

fjpotncomp = 0.0
fjpotnsurf = 0.0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

!!!    DO i = 1,ncomp
!!!        DO npt2 = 1,npot
!!!          is2 = ispot(npt2)
!!!          sum = 0.0
!!!          DO ns = 1,nsurf_sec
!!!            delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
!!!            IF (islink(ns) == is2) THEN
!!!              sum = sum - 2.0*delta_z*musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
!!!            END IF
!!!          END DO
!!!          fjpotncomp(npt2,i,jx,jy,jz) = sum     
!!!        END DO
!!!    END DO

!!!    DO is = 1,nsurf
!!!      DO npt2 = 1,npot
!!!        is2 = ispot(npt2)
!!!        sum = 0.0
!!!        IF (is == ispot(npt2)) THEN
!!!          DO ns = 1,nsurf_sec
!!!            delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
!!!            IF (islink(ns) == is2) THEN
!!!              sum = sum - 2.0*delta_z*musurf(ns,is+ncomp)*spsurf10(ns+nsurf,jx,jy,jz)
!!!            END IF
!!!          END DO
!!!          fjpotnsurf(npt2,is,jx,jy,jz) = sum
!!!        END IF   
!!!      END DO
!!!    END DO


      DO ns = 1,nsurf_sec
        surfconc = spsurf10(ns+nsurf,jx,jy,jz)
        delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

        DO i = 1,ncomp
          IF (musurf(ns,i) /= 0.0) THEN
            mutemp = musurf(ns,i)
            DO npt2 = 1,npot
             IF (ksurf(islink(ns)) == kpot(npt2)) THEN
                fjpotncomp(npt2,i,jx,jy,jz) = fjpotncomp(npt2,i,jx,jy,jz) -         &
                    2.0*delta_z*mutemp*surfconc
!!!                fjpotncomp(npt2,i,jx,jy,jz) = 0.0d0
              END IF
            END DO     
          END IF
        END DO

        DO is = 1,nsurf
          IF (musurf(ns,is+ncomp) /= 0.0) THEN
            mutemp = musurf(ns,is+ncomp)
            DO npt2 = 1,npot
              IF (ksurf(islink(ns)) == kpot(npt2)) THEN
                fjpotnsurf(npt2,is,jx,jy,jz) = fjpotnsurf(npt2,is,jx,jy,jz) -       & 
                   2.0*delta_z*mutemp*surfconc
!!!                fjpotnsurf(npt2,is,jx,jy,jz) = 0.0d0
              END IF
            END DO
          END IF   
        END DO
      END DO

    END DO
  END DO
END DO

RETURN
END SUBROUTINE jacpotential
