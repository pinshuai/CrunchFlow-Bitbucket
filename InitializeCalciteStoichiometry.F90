SUBROUTINE InitializeCalciteStoichiometry(nx,ny,nz)
USE crunchtype
USE mineral
USE concentration

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz


INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      muUranium234Bulk(jx,jy,jz) = mumin(1,kUCalcite,ik234U)
      muUranium238Bulk(jx,jy,jz) = mumin(1,kUCalcite,ik238U)
      muCalciumBulk(jx,jy,jz) = mumin(1,kUCalcite,ikCa)

    END DO
  END DO
END DO



END SUBROUTINE InitializeCalciteStoichiometry
