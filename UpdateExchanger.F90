SUBROUTINE UpdateExchanger(nx,ny,nz,nexchange)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz
INTEGER(I4B), INTENT(IN)                                   :: nexchange

!!  Internal arrays and variables

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz
INTEGER(I4B)                                               :: k

DO ix = 1,nexchange
  IF (iexchange(ix) == 1) THEN
    k = kexch(ix)
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          IF (volfx(k,jx,jy,jz) <= 0.0d0) THEN
            exchangesites(ix,jx,jy,jz) = 1.0d-30
          ELSE
            exchangesites(ix,jx,jy,jz) = cec(ix,jinit(jx,jy,jz))*wtmin(k)*volfx(k,jx,jy,jz)/(volmol(k))  
          END IF 
        END DO
      END DO
    END DO
  END IF
END DO

RETURN
END SUBROUTINE UpdateExchanger