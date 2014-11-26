SUBROUTINE surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
USE crunchtype
USE concentration
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nco


!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents
REAL(DP)                                                    :: CheckWrite

DO ns = 1,nsurf_sec

  IF (nptlink(ns) /= 0) THEN

    delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

!!  Aqueous species
    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
      sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
      sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))
    IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.0d0
    ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.693147180559945
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
      LogTotalEquivalents = LogTotalSites - 1.09861228866811
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
      LogTotalEquivalents = 0.0d0
    ELSE
      write(*,*)
      write(*,*) ' Do not recognize denticity in surface complexation'
      write(*,*)
      read(*,*)
      stop
    END IF

!!  Surface complexes
    DO is = 1,nsurf
      activity = spsurftmp(is) - LogTotalSites
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO

    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) -              &
       delta_z*2.0*LogPotential_tmp(nptlink(ns)) + sum + LogTotalEquivalents
    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))

  ELSE                                  !! Non-electrostatic option

!!  Aqueous species
    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
        sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))
    IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.0d0
    ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.693147180559945
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
      LogTotalEquivalents = LogTotalSites - 1.09861228866811
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
      LogTotalEquivalents = 0.0d0
    ELSE
      write(*,*)
      write(*,*) ' Do not recognize denticity in surface complexation'
      write(*,*)
      read(*,*)
      stop
    END IF

!!  Surface complexes
    DO is = 1,nsurf
      activity = spsurftmp(is) - LogTotalSites
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO

    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum + LogTotalEquivalents
!!    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum + LogTotalSites
    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))
  END IF

END DO

RETURN
END SUBROUTINE surf_init
!  **************************************************************
