SUBROUTINE CalciteBulkStoichiometry(nx,ny,nz,dt)
USE crunchtype
USE mineral
USE concentration

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

REAL(DP), INTENT(IN)                                       :: dt

!!  Internal arrays and variables

REAL(DP)                                                   :: volcheck
REAL(DP)                                                   :: MineralSum
REAL(DP)                                                   :: rdecay
REAL(DP)                                                   :: tot

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      volcheck = volmol(kUCalcite)*dppt(kUCalcite,jx,jy,jz)*dt
      MineralSum = volcheck + volfx(kUCalcite,jx,jy,jz)

      muUranium234Bulk(jx,jy,jz) = muUranium234Bulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum + volcheck*muUranium234(jx,jy,jz)/MineralSum
      muUranium238Bulk(jx,jy,jz) = muUranium238Bulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum + volcheck*muUranium238(jx,jy,jz)/MineralSum
      muCalciumBulk(jx,jy,jz)    = muCalciumBulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum    + volcheck*muCalcium(jx,jy,jz)/MineralSum

      rdecay = -2.88e-06*muUranium234Bulk(jx,jy,jz)*dt  + 1.55e-10*muUranium238Bulk(jx,jy,jz)*dt   !!  234U
      muUranium234Bulk(jx,jy,jz) = muUranium234Bulk(jx,jy,jz) + rdecay

      rdecay = -1.55e-10*muUranium238Bulk(jx,jy,jz)*dt                                                !!  238U
      muUranium238Bulk(jx,jy,jz) = muUranium238Bulk(jx,jy,jz) + rdecay

    END DO
  END DO
END DO


END SUBROUTINE CalciteBulkStoichiometry
