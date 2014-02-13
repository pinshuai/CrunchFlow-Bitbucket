SUBROUTINE jacsurf_init(ncomp,nsurf,nsurf_sec,nco)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nsurf
INTEGER(I4B), INTENT(IN)                             :: nsurf_sec
INTEGER(I4B), INTENT(IN)                             :: nco

!  Internal variables

REAL(DP)                                             :: sum

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: is

DO i = 1,ncomp
!        if (itype(i,nco).eq.1 .or. itype(i,nco).eq.2) then
  IF (itype(i,nco) == 1) THEN
    DO i2 = 1,ncomp+nsurf
      sum = 0.0d0
      IF (equilibrate(i,nco)) THEN
        DO ns = 1,nsurf_sec
          sum = sum + musurf(ns,i)*musurf(ns,i2)*spsurftmp10(ns+nsurf)
        END DO
      END IF
      fsurftmp(i,i2) = sum
    END DO
  ELSE
    DO i2 = 1,ncomp+nsurf
      fsurftmp(i,i2) = 0.0d0
    END DO
  END IF
END DO

DO is = 1,nsurf
  DO i2 = 1,ncomp+nsurf
    sum = 0.0
    DO ns = 1,nsurf_sec
      sum = sum + musurf(ns,is+ncomp)*musurf(ns,i2)*spsurftmp10(ns+nsurf)
    END DO
    fsurftmp(is+ncomp,i2) = sum
  END DO
  fsurftmp(is+ncomp,is+ncomp) = fsurftmp(is+ncomp,is+ncomp) + spsurftmp10(is)
END DO


RETURN
END SUBROUTINE jacsurf_init
!***********************************************************
