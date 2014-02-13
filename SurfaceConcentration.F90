SUBROUTINE SurfaceConcentration(jx,jy,jz,ncomp,nsurf,nsurf_sec)
USE crunchtype
USE concentration
USE RunTime
USE temperature
USE transport
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec

!  Internal variables and arrays

REAL(DP)                                                    :: AqueousToBulk 
INTEGER(I4B)                                                :: i,ns

!! NOTE:  Assumes a "bare" primary species site (cycling only "secondary" surface complexes)

CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

SurfaceCon = 0.0

DO i = 1,ncomp  
  DO ns = 1,nsurf_sec
    SurfaceCon(i) = SurfaceCon(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
  END DO

  IF (AqueousToBulk /= 0.0) THEN        !!  Completely dry cell
    SurfaceCon(i) = SurfaceCon(i)/AqueousToBulk
  END IF
END DO

RETURN
END SUBROUTINE SurfaceConcentration
