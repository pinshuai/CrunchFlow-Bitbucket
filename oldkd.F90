! written by S. Molins
! March 12, 2015
SUBROUTINE oldkd(ncomp,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE medium
USE transport
USE temperature

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz

!  Internal variables and arrays

INTEGER(I4B)                                  :: i
real(dp)                                      :: Retardation

DO i = 1,ncomp
  Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)
  skdold(i,jx,jy,jz) = por(jx,jy,jz) * ro(jx,jy,jz) * xgram(jx,jy,jz) * satliq(jx,jy,jz) * s(i,jx,jy,jz) * Retardation*distrib(i)   
END DO

RETURN
END SUBROUTINE
!  **************************************************************
