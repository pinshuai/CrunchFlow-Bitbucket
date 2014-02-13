!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:43:30
 
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

SUBROUTINE rocalc(tc,rotemp,nchem,DensityModule,unitsflag)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables and arrays

REAL(DP), INTENT(IN)                                  :: tc
REAL(DP), INTENT(OUT)                                 :: rotemp
INTEGER(I4B), INTENT(IN)                              :: nchem
INTEGER(I4B), DIMENSION(:), INTENT(IN)                :: unitsflag
CHARACTER (LEN=mls)                                   :: DensityModule


!  Internal variables and arrays

REAL(DP), PARAMETER                                   :: a1 = -0.004539  
REAL(DP), PARAMETER                                   :: a2 = -0.0001638  
REAL(DP), PARAMETER                                   :: a3 = 0.00002551
REAL(DP), PARAMETER                                   :: c1 = -9.9559    
REAL(DP), PARAMETER                                   :: c2 = 7.0845     
REAL(DP), PARAMETER                                   :: c3 = 3.9093
REAL(DP), PARAMETER                                   :: capa = -3.033405 
REAL(DP), PARAMETER                                   :: capb = 10.128163 
REAL(DP), PARAMETER                                   :: capc = -8.750567 
REAL(DP), PARAMETER                                   :: capd = 2.663107

REAL(DP)                                              :: term1
REAL(DP)                                              :: term2
REAL(DP)                                              :: term3
REAL(DP)                                              :: cmolal
REAL(DP)                                              :: xsum

INTEGER(I4B)                                          :: ls

IF (DensityModule == 'temperature') THEN
! Data from Helgeson and Kirkham (1974), for 500 bars
!!  rotemp = 0.102116058D+04 - 0.494459726D-01*tc -  &
!!    0.486833333D-02*tc*tc + 0.143363552D-04*tc*tc*tc - 0.219537985D-07*(tc)**04
!!rotemp = 997

!!McCutcheon, S.C., Martin, J.L, Barnwell, T.O. Jr. 1993. Water 
!!Quality in Maidment, D.R. (Editor). Handbood of Hydrology, 
!!McGraw-Hill, New York, NY (p. 11.3 )

  rotemp = 1000.0d0*(1.0d0 - (tc + 288.9414d0) / (508929.2d0*(tc + 68.12963d0))*(tc-3.9863d0)**2.0d0)


ELSE IF (DensityModule == 'sodium_chloride') THEN
  IF (unitsflag(nchem) == 5) THEN
    WRITE(*,*) 
    WRITE(*,*) ' Sodium chloride density calculation cannot be used with molarity units'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  cmolal = SQRT( ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
! Bethke's sodium chloride fit
  term1 = c1 * EXP(a1*cmolal)
  term2 = c2 * EXP(a2*tc)
  term3 = c3 * EXP(a3)          !  Assumes 1 bar pressure for now
!!  term3 = c3 * EXP(a3*pbar)   !  Original, giving dependence on pressure in bars
  xsum = term1 + term2 + term3
  rotemp = 1000.0*( capa + xsum * (capb + xsum * (capc + xsum * capd)) )

ELSE IF (DensityModule == 'sodium_nitrate') THEN
  IF (unitsflag(nchem) == 5) THEN     !  Use linear regression based on molarity from CRC
    cmolal = SQRT( ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(1.001238 + 0.051728*cmolal)
  ELSE                                !  Use linear regression based on molality from CRC
    cmolal = SQRT( ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(1.008974 + 0.04195*cmolal)
  END IF
ELSE IF (DensityModule == 'potassium_nitrate') THEN
  IF (unitsflag(nchem) == 5) THEN     !  Use linear regression based on molarity from CRC
    cmolal = SQRT( ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(0.999456119 + 0.059478795*cmolal)
  ELSE                                !  Use linear regression based on molality from CRC
    cmolal = SQRT( ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(1.002226618 + 0.052782279*cmolal)
  END IF
ELSE IF (DensityModule == 'calcium_nitrate') THEN
  IF (unitsflag(nchem) == 5) THEN     !  Use linear regression based on molarity from CRC
    cmolal = SQRT( 0.5*ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(1.001238 + 0.051728*cmolal)
  ELSE                                !  Use linear regression based on molality from CRC
    cmolal = SQRT( 0.5*ctot(MeanSalt(1),nchem)*ctot(MeanSalt(2),nchem)  )
    rotemp = 1000.0*(1.008974 + 0.04195*cmolal)
  END IF
ELSE IF (DensityModule == 'martin') THEN
  rotemp = 1000.0 + 25.0*ctot(MeanSalt(1),nchem)*1000.0
ELSE
  CALL stringlen(DensityModule,ls)
  WRITE(*,*)
  WRITE(*,*) ' Density module not recognized: ', DensityModule(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

RETURN
END SUBROUTINE rocalc
!  ***********************************************************
