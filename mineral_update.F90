!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:52:46
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE mineral_update(nx,ny,nz,nrct,dt,dtnew,ineg,jpor,deltmin)
USE crunchtype
USE params
USE runtime, ONLY: voltol,inagaki,inagaki2,ncounter,ReadGautier
USE concentration
USE mineral
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: nrct
INTEGER(I4B), INTENT(IN)                                    :: jpor
INTEGER(I4B), INTENT(OUT)                                   :: ineg

REAL(DP), INTENT(IN)                                        :: dt
REAL(DP), INTENT(INOUT)                                     :: dtnew
REAL(DP), INTENT(IN)                                        :: deltmin

!  Internal variables and arrays

REAL(DP)                                                    :: volcheck
REAL(DP)                                                    :: voldiff
REAL(DP)                                                    :: volneg
REAL(DP)                                                    :: vinit
REAL(DP)                                                    :: sum
REAL(DP)                                                    :: porfactor
REAL(DP)                                                    :: volumeterm
REAL(DP)                                                    :: dtTemp

INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
INTEGER(I4B)                                                :: jxneg
INTEGER(I4B)                                                :: jyneg
INTEGER(I4B)                                                :: jzneg
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: ireduce
INTEGER(I4B)                                                :: kchg
INTEGER(I4B)                                                :: kneg
INTEGER(I4B)                                                :: ls

!!  Hardwired for ripening/densification reaction
INTEGER(I4B)                                                :: kJapaneseGlassSilicaOnly
INTEGER(I4B)                                                :: kJapaneseGlassNoSilica
INTEGER(I4B)                                                :: kJapaneseGlassSilicaGel
INTEGER(I4B)                                                :: kSiO2Amorphous
INTEGER(I4B)                                                :: kSiO2Dense
INTEGER(I4B)                                                :: iTracerH2O
REAL(DP)                                                    :: RateDensify 
REAL(DP)                                                    :: inhibition
REAL(DP)                                                    :: catalysis
REAL(DP)                                                    :: rdensify
REAL(DP)                                                    :: DensificationFactor

CHARACTER (LEN=mls)                                         :: dumstring

INTEGER(I4B)                                                :: nn 
INTEGER(I4B)                                                :: np


!!  ********* Hardwired for ripening/densification reaction
kJapaneseGlassSilicaOnly = 1
kJapaneseGlassNoSilica = 2  
kJapaneseGlassSilicaGel = 16
kSiO2Dense = 3
kSiO2Amorphous = 2
iTracerH2O = 16
!!RateDensify  = 750.0d0
RateDensify = 0.0d0

!!  ********** Hardwired for ripening/densification reaction

ineg = 0
ireduce = 0
volneg = 0.0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO k = 1,nrct
!!       IF (volsave(k,jx,jy,jz) == 0.0d0) then
!!          volsave(k,jx,jy,jz) = volfx(k,jx,jy,jz)
!!        END IF
        IF (MineralAssociate(k)) THEN
          kk = MineralID(k)
          volcheck = volmol(kk)*dppt(k,jx,jy,jz)*dt   !  Point to volume fraction of associated mineral
          voldiff = volcheck + volfx(kk,jx,jy,jz)
          IF (voldiff < volneg) THEN
            volneg = voldiff
            kneg = kk
            jxneg = jx
            jyneg = jy
            jzneg = jz
          END IF
        ELSE
          volcheck = volmol(k)*dppt(k,jx,jy,jz)*dt
          voldiff = volcheck + volfx(k,jx,jy,jz)
          IF (voldiff < volneg) THEN
            volneg = voldiff
            kneg = k
            jxneg = jx
            jyneg = jy
            jzneg = jz
          END IF
        END IF


      END DO
    END DO
  END DO
END DO

!300  IF (ireduce == 1) THEN
!  dtnew = dt/2.0
!  dumstring = umin(kchg)
!  call stringlen(dumstring,ls)
!  WRITE(*,*) ' Mineral dissolving too fast: ', dumstring(1:ls)
!  WRITE(*,100) dtnew
!  WRITE(*,*)
!  ineg = 1
!  RETURN
!END IF
!IF (ireduce == 2) THEN
!  dtnew = dt/2.0
!  dumstring = umin(kchg)
!  call stringlen(dumstring,ls)
!  WRITE(*,*) ' Mineral precipitating too fast: ', dumstring(1:ls)
!  WRITE(*,100) dtnew
!  WRITE(*,*)
!  ineg = 1
!  RETURN
!END IF

!  Calculate a new delt and repeat time step if mineral concentrations
!  are too far below 0

IF (volneg < -voltol) THEN
  dtTemp = volfx(kneg,jxneg,jyneg,jzneg)/(-dppt(kneg,jxneg,jyneg,jzneg)*volmol(kneg))
  IF (dtTemp > deltmin) THEN
    dtnew = dtTemp
  ELSE
    dtnew = deltmin
  END IF
  dumstring = umin(kneg)
  call stringlen(dumstring,ls)
  WRITE(*,*) ' Mineral volume fraction less than zero: ', dumstring(1:ls)
  WRITE(*,*) ' New time step: ',dtnew
  WRITE(*,*)
  ineg = 1
  RETURN
END IF

ncounter = ncounter + 1

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        
      DO k = 1,nrct

        IF (MineralAssociate(k)) THEN
            
          kk = MineralID(k)
          volcheck = volmol(kk)*dppt(k,jx,jy,jz)*dt   !  Point to volume fraction of associated mineral
          volfx(kk,jx,jy,jz) = volfx(kk,jx,jy,jz) + volcheck
          IF (volfx(kk,jx,jy,jz) < 0.0) THEN
            volfx(kk,jx,jy,jz) = 0.0
          END IF
          
        ELSE
            
          volcheck = volmol(k)*dppt(k,jx,jy,jz)*dt
          volfx(k,jx,jy,jz) = volfx(k,jx,jy,jz) + volcheck
          volSaveByTimeStep(ncounter,k,jx,jy,jz) = volcheck
          volSave(k,jx,jy,jz) = volSave(k,jx,jy,jz) + volSaveByTimeStep(ncounter,k,jx,jy,jz)

          IF (ncounter > 100) then
            volSave(k,jx,jy,jz) = volSave(k,jx,jy,jz) - volSaveByTimeStep(1,k,jx,jy,jz)
            do nn = ncounter,2,-1
              volSaveByTimeStep(nn-1,k,jx,jy,jz) = volSaveByTimeStep(nn,k,jx,jy,jz)
            end do
            ncounter = ncounter - 1
          END IF

          IF (volfx(k,jx,jy,jz) < 0.0) THEN
            volfx(k,jx,jy,jz) = 0.0
          END IF

        END IF

      END DO
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      IF (inagaki) THEN

!!  ********* Hardwired for ripening/densification reaction
       inhibition =  (0.01d0/(volfx(kJapaneseGlassNoSilica,jx,jy,jz) + 0.01d0))**2.0
!!       catalysis  = volfx(kJapaneseGlassSilicaOnly,jx,jy,jz)/                 &
!!            (volfx(kJapaneseGlassSilicaOnly,jx,jy,jz) + 0.01)
       catalysis = volfx(kJapaneseGlassSilicaOnly,jx,jy,jz)
       rdensify = RateDensify*s(iTracerH2O,jx,jy,jz)*catalysis*inhibition
       volfx(kJapaneseGlassSilicaGel,jx,jy,jz) = volfx(kJapaneseGlassSilicaGel,jx,jy,jz) &
                  + 1.50*rdensify*dt 
!!                  + 2.105d0*rdensify*dt
       volfx(kJapaneseGlassSilicaOnly,jx,jy,jz) = volfx(kJapaneseGlassSilicaOnly,jx,jy,jz) &
                  - rdensify*dt


!!  ********** Hardwired for ripening/densification reaction

      END IF

      IF (inagaki2) THEN

!!         write(*,*) ' Densifying....'
!!         read(*,*)
!!  ********* Hardwired for ripening/densification reaction
!!       inhibition =  (0.01d0/(volfx(kJapaneseGlassNoSilica,jx,jy,jz) + 0.01d0))**2.0
!!       catalysis  = volfx(kJapaneseGlassSilicaOnly,jx,jy,jz)/                 &
!!            (volfx(kJapaneseGlassSilicaOnly,jx,jy,jz) + 0.01)

       catalysis = volfx(kSiO2Amorphous,jx,jy,jz)
       inhibition = 1.0d0
       rdensify = RateDensify*s(iTracerH2O,jx,jy,jz)*catalysis*inhibition
       DensificationFactor = volmol(kSiO2Dense)/volmol(kSiO2Amorphous)
       volfx(kSiO2Dense,jx,jy,jz) = volfx(kSiO2Dense,jx,jy,jz) &
                  + DensificationFactor*rdensify*dt 
       volfx(kSiO2Amorphous,jx,jy,jz) = volfx(kSiO2Amorphous,jx,jy,jz) &
                  - rdensify*dt


!!  ********** Hardwired for ripening/densification reaction

      END IF

      sum = 0.0
      DO k = 1,nrct

        vinit = volin(k,jinit(jx,jy,jz))
        
        IF (LocalEquilibrium(k)) THEN                      !! Local equilibrium fantasy, so don't change the surface area
          IF (volfx(k,jx,jy,jz) > 0.0d0) THEN
            area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))
          ELSE
            area(k,jx,jy,jz) = 0.0d0
          END IF
          sum = sum + volfx(k,jx,jy,jz)
        ELSE                                               !! Update reactive surface area

          IF (iarea(k,jinit(jx,jy,jz)) == 0) THEN                     !!  Bulk surface area
              
            IF (vinit == 0.0d0) THEN
              area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))* (volfx(k,jx,jy,jz)/0.01)**0.6666
              sum = sum + volfx(k,jx,jy,jz)
            ELSE
              area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))* (volfx(k,jx,jy,jz)/vinit)**0.6666
              sum = sum + volfx(k,jx,jy,jz)
            END IF
            
          ELSE                                                        !!  Specific surface area
              
            IF ( volin(k,jinit(jx,jy,jz)) == 0.0d0 .AND. volfx(k,jx,jy,jz) < voltemp(k,jinit(jx,jy,jz)) ) THEN   !!  Initially a zero volume fraction, so use "threshold" volume fraction
              area(k,jx,jy,jz) = voltemp(k,jinit(jx,jy,jz))*specific(k,jinit(jx,jy,jz))*wtmin(k)/volmol(k)
            ELSE
              area(k,jx,jy,jz) = volfx(k,jx,jy,jz)*specific(k,jinit(jx,jy,jz))*wtmin(k)/volmol(k)
            END IF
            sum = sum + volfx(k,jx,jy,jz)
            
          END IF
          
          DO np = 1,nreactmin(k)
            IF (imintype(np,k) == 10) THEN
              area(k,jx,jy,jz) = volfx(k,jx,jy,jz)*SurfaceAreaNucleation(np,k)*wtmin(k)/volmol(k)
              IF (volfx(k,jx,jy,jz) > 0.0d0) THEN
                CONTINUE
              END IF
            END IF
          END DO

        END IF

      END DO
      
!!  Update porosity, with a save of the porosity to porOld

      porold(jx,jy,jz) = por(jx,jy,jz)
      IF (jpor /= 1 .AND. .NOT. ReadGautier) THEN
        por(jx,jy,jz) = porin(jx,jy,jz)
      ELSE
        por(jx,jy,jz) = 1.0 - sum
        IF (por(jx,jy,jz) < MinimumPorosity) THEN
          por(jx,jy,jz) = MinimumPorosity
        END IF
      END IF
      
      
      porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**0.6666666666666
!!      IF (porfactor < 1.0d0) THEN
!!        DO k = 1,nrct
!!          area(k,jx,jy,jz) = area(k,jx,jy,jz)*porfactor
!!        END DO
!!      END IF
    END DO
  END DO
END DO



200 FORMAT(2X,'Porosity has gone to zero')
201 FORMAT(2X,a18,2X,1PE12.4,2X,'At grid point ',i2)
100 FORMAT(1X,'Time step reduced to ',1PE12.4)

RETURN
END SUBROUTINE mineral_update
!  *******************************************************
