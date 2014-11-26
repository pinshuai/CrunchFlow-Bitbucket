SUBROUTINE oldsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
!USE concentration, ONLY: ssurfn,ksurf,jinit
!USE mineral, ONLY: volfx,volmol,wtmin,site_density,specific
USE concentration
USE mineral
USE medium
USE temperature

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz

!  Internal variables and arrays

REAL(DP)                                      :: sum
REAL(DP)                                      :: permole
REAL(DP)                                      :: check
REAL(DP)                                                    :: volMinimum

INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: k
INTEGER(I4B)                                  :: ns

DO is = 1,nsurf

  k = ksurf(is)
  permole = site_density(is,jinit(jx,jy,jz))*specific(k,jinit(jx,jy,jz))*wtmin(k)    !  Mole sites/Mole mineral

  IF (volin(k,jinit(jx,jy,jz)) == 0.0d0 .AND. volfx(k,jx,jy,jz)< voltemp(k,jinit(jx,jy,jz)) ) THEN
    ssurfn(is,jx,jy,jz) = permole*voltemp(k,jinit(jx,jy,jz))/(volmol(k))
  ELSE
    volMinimum = volfx(k,jx,jy,jz)
    if (volMinimum < 1.0D-15) then
       volMinimum = 1.0D-15
    end if
    ssurfn(is,jx,jy,jz) = permole*volMinimum/(volmol(k))
  END IF

!!  IF (ssurfn(is,jx,jy,jz) < 1.D-20) THEN
!!    ssurfn(is,jx,jy,jz) = 1.D-20
!!  END IF

  LogTotalSurface(is,jx,jy,jz) = DLOG(ssurfn(is,jx,jy,jz))

END DO

RETURN
END SUBROUTINE oldsurf
!  **************************************************************
