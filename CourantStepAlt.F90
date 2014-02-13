SUBROUTINE CourantStepAlt(nx,ny,nz,dtmaxcour)
USE crunchtype
USE runtime
USE medium
USE transport
USE flow
USE temperature
USE modflowModule
USE ReadFlow

IMPLICIT NONE    

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

REAL(DP), INTENT(OUT)                                    :: dtmaxcour

!  Internal variables and arrays

INTEGER(I4B)                                             :: jx
INTEGER(I4B)                                             :: jy
INTEGER(I4B)                                             :: jz
INTEGER(I4B)                                             :: i

REAL(DP)                                                 :: dtemp
REAL(DP)                                                 :: porsatro
REAL(DP)                                                 :: dtmin

IF (xflow .OR. yflow .OR. zflow) THEN

  IF (ReadNuft) THEN        !  Mass flux (kg/m**2/yr) read from NUFT
    dtmin = 1000000.0

    DO jz = 0,nz
      DO jy = 0,ny
        DO jx = 0,nx
          IF (jy /= 0 .AND. jz /= 0) THEN
            IF (qx(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dxx(jx)*porsatro/ABS(qx(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jz /= 0) THEN
            IF (qy(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dyy(jy)*porsatro/ABS(qy(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jy /= 0) THEN
            IF (qz(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dzz(jx,jy,jz)*porsatro/ABS(qz(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
        END DO
      END DO
    END DO

!   Check the recharge for the Courant criteria

    DO jy=1,ny
      DO jx=1,nx
        porsatro = por(jx,jy,1)*satliq(jx,jy,1)*ro(jx,jy,1)
        IF(qrecharge(jx,jy) == 0.0)THEN
          dtemp=1.e30
        ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
          dtemp=ABS( 0.5*dxx(jx)*dyy(jy)*dzz(jx,jy,1)*porsatro/qrecharge(jx,jy)  )
!fp! end_onproc();
        END IF
        dtmin=MIN(dtemp,dtmin)
      END DO
    END DO

    DO i = 1, NumSourceTerms
      jx = jxNuftSource(i)
      jy = jyNuftSource(i)
      jz = jzNuftSource(i)
      porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
      IF(qg(jx,jy,jz) == 0.0)THEN
        dtemp=1.e30
      ELSE
!fp! if_onproc( {#expr# por(jx,jy,jz) #});
        dtemp=ABS(0.5*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*porsatro/qg(jx,jy,jz))
!fp! end_onproc();
      END IF
      dtmin=MIN(dtemp,dtmin)
    END DO

    dtmaxcour = dtmin

!fp! reduce( {#ident# dtmaxcour #}, min_op);

  ELSE                   !  Volumetric flux used

    dtmin = 1000000.0
    DO jz = 0,nz
      DO jy = 0,ny
        DO jx = 0,nx
          IF (jy /= 0 .AND. jz /= 0) THEN
            IF (qx(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dxx(jx)*porsatro/ABS(qx(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jz /= 0) THEN
            IF (qy(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dyy(jy)*porsatro/ABS(qy(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jy /= 0) THEN
            IF (qz(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dzz(jx,jy,jz)*porsatro/ABS(qz(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
        END DO
      END DO
    END DO
    dtmaxcour = dtmin

!fp! reduce( {#ident# dtmaxcour #}, min_op);

  END IF
END IF 

RETURN
END SUBROUTINE CourantStepAlt
