subroutine FluxWeightedConcentrationWrite(ncomp,nspec,nx,ny,nz,nn,time,delt )
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE strings
USE ReadFlow
USE modflowModule
USE NanoCrystal

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz
INTEGER(I4B), INTENT(IN)                                     :: nn

REAL(DP), INTENT(IN)                                         :: time
REAL(DP), INTENT(IN)                                         :: delt

!! INTERNAL VARIABLES

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: is
INTEGER(I4B)                                                  :: nex
INTEGER(I4B)                                                  :: j
INTEGER(I4B)                                                  :: nxyz
INTEGER(I4B)                                                  :: jx
INTEGER(I4B)                                                  :: jy
INTEGER(I4B)                                                  :: jz
INTEGER(I4B)                                                  :: ll
INTEGER(I4B)                                                  :: intfile
INTEGER(I4B)                                                  :: jx_lo 
INTEGER(I4B)                                                  :: jx_hi
INTEGER(I4B)                                                  :: jy_lo 
INTEGER(I4B)                                                  :: jy_hi
INTEGER(I4B)                                                  :: jz_lo 
INTEGER(I4B)                                                  :: jz_hi

CHARACTER (LEN=mls)                                           :: filename
CHARACTER (LEN=mls)                                           :: dummy
CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=12)                                            :: writeph

LOGICAL(LGT)                                                  :: ext

LOGICAL(LGT)                                                  :: X_CapturePlane
LOGICAL(LGT)                                                  :: y_CapturePlane
LOGICAL(LGT)                                                  :: Z_CapturePlane

REAL(DP)                                                      :: eps
REAL(DP)                                                      :: PrintTime
REAL(DP)                                                      :: qxSum
REAL(DP)                                                      :: xflux
REAL(DP)                                                      :: qySum
REAL(DP)                                                      :: yflux
REAL(DP)                                                      :: qzSum
REAL(DP)                                                      :: zflux



!!  Write to file at fixed intervals
  
IF (nplotFluxWeightedConcentration > 0) THEN

  IF (MOD(nn,interval) == 0) THEN

    PrintTime = OutputTimeScale*time

    DO ll = 1,nFluxWeightedConcentrationFile

      X_CapturePlane = .FALSE.
      Y_CapturePlane = .FALSE.
      Z_CapturePlane = .FALSE.

      intfile = 75+ll
      jx_lo = jxFluxWeightedConcentration_lo(ll)
      jx_hi = jxFluxWeightedConcentration_hi(ll)
      jy_lo = jyFluxWeightedConcentration_lo(ll)
      jy_hi = jyFluxWeightedConcentration_hi(ll)
      jz_lo = jzFluxWeightedConcentration_lo(ll)
      jz_hi = jzFluxWeightedConcentration_hi(ll)

!!   Check to see where this is an X, Y, or Z capture plane

      IF (nx > 1 .AND. jx_lo == jx_hi) THEN
        X_CapturePlane = .TRUE.
        Y_CapturePlane = .FALSE.
        Z_CapturePlane = .FALSE.
      END IF

      IF (ny > 1 .AND. jy_lo == jy_hi .AND. X_CapturePlane==.FALSE.) THEN
        Y_CapturePlane = .TRUE.
        Z_CapturePlane = .FALSE.
      END IF

      IF (nz > 1 .AND. jz_lo == jz_hi .AND. X_CapturePlane==.FALSE. .AND. Y_CapturePlane==.FALSE.) THEN
        Z_CapturePlane = .TRUE.
      END IF


      IF (tecplot) THEN

!!    NOTE:  Based on the Right Hand Rule, with fluxes evaluated at the Right Hand Side

        IF (X_CapturePlane) THEN

          qxSum = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi
                qxSum = qxSum + qx(jx,jy,jz)
              END DO
            END DO
          END DO

          XfluxWeightedConcentration = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi

                DO i = 1,nplotFluxWeightedConcentration
                  IF (qx(jx,jy,jz) > 0.0d0) THEN
                    xflux = (qx(jx,jy,jz)/qxSum)*s(iplotFluxWeightedConcentration(i),jx,jy,jz)
                  ELSE
                    xflux = (qx(jx,jy,jz)/qxSum)*s(iplotFluxWeightedConcentration(i),jx+1,jy,jz)
                  END IF
                  XfluxWeightedConcentration(i) = XfluxWeightedConcentration(i) + xflux
                END DO

              END DO
            END DO
          END DO
              
          WRITE(intfile,185) PrintTime,(XfluxWeightedConcentration(i),i=1,nplotFluxWeightedConcentration)     

        ELSE IF (Y_CapturePlane) THEN

          qySum = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi
                qySum = qySum + qy(jx,jy,jz)
              END DO
            END DO
          END DO

          YfluxWeightedConcentration = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi

                DO i = 1,nplotFluxWeightedConcentration
                  IF (qy(jx,jy,jz) > 0.0d0) THEN
                    yflux = (qy(jx,jy,jz)/qySum)*s(iplotFluxWeightedConcentration(i),jx,jy,jz)
                  ELSE
                    yflux = (qy(jx,jy,jz)/qySum)*s(iplotFluxWeightedConcentration(i),jx,jy+1,jz)
                  END IF
                  YfluxWeightedConcentration(i) = YfluxWeightedConcentration(i) + yflux
                END DO

              END DO
            END DO
          END DO
              
          WRITE(intfile,185) PrintTime,(YfluxWeightedConcentration(i),i=1,nplotFluxWeightedConcentration)     

        ELSE IF (Z_CapturePlane) THEN

          qzSum = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi
                qzSum = qzSum + qz(jx,jy,jz)
              END DO
            END DO
          END DO

          ZfluxWeightedConcentration = 0.0d0
          DO jz = jz_lo,jz_hi
            DO jy = jy_lo,jy_hi
              DO jx = jx_lo,jx_hi

                DO i = 1,nplotFluxWeightedConcentration
                  IF (qz(jx,jy,jz) > 0.0d0) THEN
                    zflux = (qz(jx,jy,jz)/qzSum)*s(iplotFluxWeightedConcentration(i),jx,jy,jz)
                  ELSE
                    zflux = (qz(jx,jy,jz)/qzSum)*s(iplotFluxWeightedConcentration(i),jx,jy,jz+1)
                  END IF
                  ZfluxWeightedConcentration(i) = ZfluxWeightedConcentration(i) + zflux
                END DO

              END DO
            END DO
          END DO
              
          WRITE(intfile,185) PrintTime,(ZfluxWeightedConcentration(i),i=1,nplotFluxWeightedConcentration)    

        END IF       
                         
      ELSE    !! Kaleidagraph (and OriginLab?)

!!    NOTE:  Based on the Right Hand Rule, with fluxes evaluated at the Right Hand Side

        qxSum = 0.0d0
        DO jz = jz_lo,jz_hi
          DO jy = jy_lo,jy_hi
            DO jx = jx_lo,jx_hi
              qxSum = qxSum + qx(jx,jy,jz)
            END DO
          END DO
        END DO

        XfluxWeightedConcentration = 0.0d0
        DO jz = jz_lo,jz_hi
          DO jy = jy_lo,jy_hi
            DO jx = jx_lo,jx_hi

              DO i = 1,nplotFluxWeightedConcentration
                IF (qx(jx,jy,jz) > 0.0d0) THEN
                  xflux = (qx(jx,jy,jz)/qxSum)*s(iplotFluxWeightedConcentration(i),jx,jy,jz)
                ELSE
                  xflux = (qx(jx,jy,jz)/qxSum)*s(iplotFluxWeightedConcentration(i),jx+1,jy,jz)
                END IF
                XfluxWeightedConcentration(i) = XfluxWeightedConcentration(i) + xflux
              END DO

            END DO
          END DO
        END DO

       WRITE(intfile,185) PrintTime,(XfluxWeightedConcentration(i),i=1,nplotFluxWeightedConcentration)  

      END IF

    END DO

  END IF
END IF

705 FORMAT(1X,1PE12.5,13x,150(1X,1PE22.14))
185 FORMAT(1X,1PE12.5,13x,100(1X,1PE22.14))



RETURN

END subroutine FluxWeightedConcentrationWrite