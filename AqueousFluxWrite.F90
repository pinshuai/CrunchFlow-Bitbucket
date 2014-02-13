subroutine AqueousFluxWrite(ncomp,nspec,nx,ny,nz,nn,time,delt )
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

REAL(DP)                                                      :: eps
REAL(DP)                                                      :: PrintTime
REAL(DP)                                                      :: xDiffuse
REAL(DP)                                                      :: xflux

!!REAL(DP), DIMENSION(ncomp)                                    :: CumulativeXFlux

!!  Write to file at fixed intervals
  
IF (nplotAqueousFlux > 0) THEN

  IF (MOD(nn,interval) == 0) THEN

    PrintTime = OutputTimeScale*time

    DO ll = 1,nAqueousFluxSeriesFile

      intfile = 50+ll
      jx_lo = jxAqueousFluxSeries_lo(ll)
      jx_hi = jxAqueousFluxSeries_hi(ll)
      jy_lo = jyAqueousFluxSeries_lo(ll)
      jy_hi = jyAqueousFluxSeries_hi(ll)
      jz_lo = jzAqueousFluxSeries_lo(ll)
      jz_hi = jzAqueousFluxSeries_hi(ll)

      IF (tecplot) THEN

!!    NOTE:  Based on the Right Hand Rule, with fluxes evaluated at the Right Hand Side

        InstantaneousXflux = 0.0d0
        DO jz = jz_lo,jz_hi
          DO jy = jy_lo,jy_hi
            DO jx = jx_lo,jx_hi

              DO i = 1,nplotAqueousFlux
                xDiffuse = netDiffuseX(jx,jy,jz)*( s(iplotAqueousFlux(i),jx+1,jy,jz) - s(iplotAqueousFlux(i),jx,jy,jz) )
                IF (netflowX(jx,jy,jz) > 0.0d0) THEN
                  xflux = netflowX(jx,jy,jz)*s(iplotAqueousFlux(i),jx,jy,jz) + xDiffuse
                ELSE
                  xflux = netflowX(jx,jy,jz)*s(iplotAqueousFlux(i),jx+1,jy,jz) + xDiffuse
                END IF
                InstantaneousXflux(i) = InstantaneousXflux(i) + xflux
                CumulativeXflux(i,ll) = CumulativeXflux(i,ll) + xflux*delt
              END DO

            END DO
          END DO
        END DO
              
        WRITE(intfile,185) PrintTime,(CumulativeXflux(i,ll),i=1,nplotAqueousFlux)            
!!        WRITE(intfile,185) PrintTime,(CumulativeXflux(i,ll),i=1,nplotAqueousFlux),       &
!!                             (InstantaneousXflux(i),i=1,nplotAqueousFlux)
                            


      ELSE    !! Kaleidagraph (and OriginLab?)

!!    NOTE:  Based on the Right Hand Rule, with fluxes evaluated at the Right Hand Side

        InstantaneousXflux = 0.0d0
         DO jz = jz_lo,jz_hi
          DO jy = jy_lo,jy_hi
            DO jx = jx_lo,jx_hi

              DO i = 1,nplotAqueousFlux
                xDiffuse = netDiffuseX(jx,jy,jz)*( s(iplotAqueousFlux(i),jx+1,jy,jz) - s(iplotAqueousFlux(i),jx,jy,jz) )
                IF (netflowX(jx,jy,jz) > 0.0d0) THEN
                  xflux = netflowX(jx,jy,jz)*s(iplotAqueousFlux(i),jx,jy,jz) + xDiffuse
                ELSE
                  xflux = netflowX(jx,jy,jz)*s(iplotAqueousFlux(i),jx+1,jy,jz) + xDiffuse
                END IF
                InstantaneousXflux(i) = InstantaneousXflux(i) + xflux
                CumulativeXflux(i,ll) = CumulativeXflux(i,ll) + xflux*delt
              END DO

            END DO
          END DO
        END DO

        WRITE(intfile,185) PrintTime,(CumulativeXflux(i,ll),i=1,nplotAqueousFlux)             
!!        WRITE(intfile,185) PrintTime,(CumulativeXflux(i,ll),i=1,nplotAqueousFlux),       &
!!                             (InstantaneousXflux(i),i=1,nplotAqueousFlux)

      END IF

    END DO

  END IF
END IF

705 FORMAT(1X,1PE12.5,13x,150(1X,1PE22.14))
185 FORMAT(1X,1PE12.5,13x,100(1X,1PE22.14))



RETURN

END subroutine AqueousFluxWrite