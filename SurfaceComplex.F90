SUBROUTINE SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)
USE crunchtype
USE concentration
USE mineral
USE medium
USE transport
USE temperature
USE RunTime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: AqueousToBulk
REAL(DP)                                                    :: LogAqueousToBulk
REAL(DP)                                                    :: MeanSaltConcentration
REAL(DP)                                                    :: MassFraction
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

!!      CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
!!      LogAqueousToBulk = DLOG(AqueousToBulk)

      DO ns = 1,nsurf_sec

        IF (nptlink(ns) /= 0) THEN                            !  Electrostatic correction

          delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF

          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 
          IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
            LogTotalEquivalents = LogTotalSites - 0.0d0
          ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
            LogTotalEquivalents = LogTotalSites - 0.693147180559945
          ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
            LogTotalEquivalents = LogTotalSites - 1.09861228866811
          ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
            LogTotalEquivalents = 0.0d0
          ELSE
            write(*,*)
            write(*,*) ' Do not recognize denticity in surface complexation'
            write(*,*)
            read(*,*)
            stop
          END IF

          DO is = 1,nsurf

!!            activity = spsurf(is,jx,jy,jz) - LogTotalSites - LogAqueousToBulk
            activity = spsurf(is,jx,jy,jz) - LogTotalSites
            sum = sum + musurf(ns,is+ncomp)*activity

          END DO

!! NOTE:  Below is the LOg concentration of sites in units of mol/kgw
          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum -                     &
             delta_z*2.0*LogPotential(nptlink(ns),jx,jy,jz) + LogTotalEquivalents

!!          spsurf10(ns+nsurf,jx,jy,jz) = AqueousToBulk*DEXP( spsurf(ns+nsurf,jx,jy,jz) )
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        ELSE                                                  !  Non-electrostatic 

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF

          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 
          IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
            LogTotalEquivalents = LogTotalSites - 0.0d0
          ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
            LogTotalEquivalents = LogTotalSites - 0.693147180559945
          ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
            LogTotalEquivalents = LogTotalSites - 1.09861228866811
          ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
            LogTotalEquivalents = 0.0d0
          ELSE
            write(*,*)
            write(*,*) ' Do not recognize denticity in surface complexation'
            write(*,*)
            read(*,*)
            stop
          END IF

          DO is = 1,nsurf

            activity = spsurf(is,jx,jy,jz) - LogTotalSites
            sum = sum + musurf(ns,is+ncomp)*activity

          END DO

          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum + LogTotalEquivalents

!!          spsurf10(ns+nsurf,jx,jy,jz) = AqueousToBulk*DEXP( spsurf(ns+nsurf,jx,jy,jz) )
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        END IF
      END DO

    END DO
  END DO
END DO

RETURN
END SUBROUTINE SurfaceComplex
!  **************************************************************
