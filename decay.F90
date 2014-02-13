SUBROUTINE decay(ncomp,ikin,delt,jx,jy,jz)
USE crunchtype
USE params
USE mineral
USE concentration

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: ncomp
INTEGER(I4B), INTENT(IN)                               :: ikin
REAL(DP), INTENT(IN)                                   :: delt
INTEGER(I4B), INTENT(IN)                               :: jx
INTEGER(I4B), INTENT(IN)                               :: jy
INTEGER(I4B), INTENT(IN)                               :: jz

!  Internal variables and arrays

INTEGER(I4B)                                           :: ir
INTEGER(I4B)                                           :: kd
INTEGER(I4B)                                           :: k
INTEGER(I4B)                                           :: np
INTEGER(I4B)                                           :: i

REAL(DP)                                               :: rdecay
REAL(DP)                                               :: voltemporary
REAL(DP)                                               :: frx
REAL(DP)                                               :: denom

DO ir = 1,ikin

  IF (iaqtype(ir) == 4) then
    DO kd = 1,nrad_decay(ir)
      k = kradpoint(kd,ir)
      np = npradpoint(kd,ir)
      IF (iraddecay(ir) == 0) THEN
        WRITE(*,*) ' Error in decay reactions:  parent species in mineral not found'
        READ(*,*)
        STOP
      END IF
      IF (dppt(k,jx,jy,jz) > 0.0) THEN
        voltemporary = volmol(k)*dppt(k,jx,jy,jz)*delt
        denom = voltemporary + volfx(k,jx,jy,jz)
        IF (denom == 0.0) then      !  Only happens with zero molar volume
          frx = 1.0
        ELSE
          frx = voltemporary/denom
        END IF
        DO i = 1,ncomp
          mumin_decay(np,k,i,jx,jy,jz) = frx*mumin(np,k,i) + (1.0-frx)*mumin_decay(np,k,i,jx,jy,jz)
          rdecay = mukin(ir,i)*ratek(1,ir)*mumin_decay(np,k,iraddecay(ir),jx,jy,jz)*delt
          mumin_decay(np,k,i,jx,jy,jz) = mumin_decay(np,k,i,jx,jy,jz) +  rdecay
        END DO
      ELSE
        IF (volfx(k,jx,jy,jz) > 0.0) THEN
          DO i = 1,ncomp
            rdecay = mukin(ir,i)*ratek(1,ir)*mumin_decay(np,k,iraddecay(ir),jx,jy,jz)*delt
            mumin_decay(np,k,i,jx,jy,jz) = mumin_decay(np,k,i,jx,jy,jz) + rdecay
          END DO
        END IF
      END IF
    END DO
  END IF

END DO
    

END SUBROUTINE decay
