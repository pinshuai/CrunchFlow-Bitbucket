SUBROUTINE JacPotentialLocal(ncomp,nsurf,nsurf_sec,npot,jx,jy,jz)
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
INTEGER(I4B), INTENT(IN)                                                :: jx
INTEGER(I4B), INTENT(IN)                                                :: jy
INTEGER(I4B), INTENT(IN)                                                :: jz

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

fjpotncomp_local = 0.0
fjpotnsurf_local = 0.0

DO ns = 1,nsurf_sec
  surfconc = spsurf10(ns+nsurf,jx,jy,jz)
  delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

  DO i = 1,ncomp
    IF (musurf(ns,i) /= 0.0) THEN
      mutemp = musurf(ns,i)
      DO npt2 = 1,npot
        is2 = ispot(npt2)
        IF (islink(ns) == is2) THEN
          fjpotncomp_local(npt2,i) = fjpotncomp_local(npt2,i) -         &
                    2.0*delta_z*mutemp*surfconc
        END IF
      END DO     
    END IF
  END DO

  DO is = 1,nsurf
    IF (musurf(ns,is+ncomp) /= 0.0) THEN
      mutemp = musurf(ns,is+ncomp)
      DO npt2 = 1,npot
        is2 = ispot(npt2)
        IF (is == is2 .AND. islink(ns) == is2) THEN
          fjpotnsurf_local(npt2,is) = fjpotnsurf_local(npt2,is) -       & 
                   2.0*delta_z*mutemp*surfconc
        END IF
      END DO
    END IF   
  END DO

END DO


RETURN
END SUBROUTINE JacPotentialLocal
