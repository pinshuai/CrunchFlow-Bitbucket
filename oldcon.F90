SUBROUTINE oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz

!  Internal variables and arrays

REAL(DP)                                      :: sum1
REAL(DP)                                      :: sum2
REAL(DP)                                      :: sum3

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: ksp
INTEGER(I4B)                                  :: nex
INTEGER(I4B)                                  :: ns

DO i = 1,ncomp
  sum1 = 0.0
  sum2 = 0.0
  sum3 = 0.0
  DO ksp = 1,nspec
    sum1 = sum1 + muaq(ksp,i)*sp10(ksp+ncomp,jx,jy,jz)
  END DO
  DO nex = 1,nexch_sec
    sum2 = sum2 + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)   
  END DO
  DO ns = 1,nsurf_sec
    sum3 = sum3 + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
  END DO
  sn(i,jx,jy,jz) = sum1 + sp10(i,jx,jy,jz)
  sexold(i,jx,jy,jz) = sum2
  ssurfold(i,jx,jy,jz) = sum3
END DO

RETURN
END SUBROUTINE oldcon
!  **************************************************************
