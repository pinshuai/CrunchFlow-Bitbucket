!************** (C) COPYRIGHT 1995 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:55:57
 
!                      All Rights Reserved

!  OSRT (Operator Splitting Reactive Transport) IS PROVIDED "AS IS"
!  AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED. THE USER ASSUMES ALL RISKS
!  OF USING OSRT. THERE IS NO CLAIM OF THE MERCHANTABILITY OR FITNESS
!  FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

SUBROUTINE density(jx,jy,jz)
USE crunchtype
USE params
USE temperature
USE concentration
USE RunTime

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                :: jz
INTEGER(I4B), INTENT(IN)                                :: jx
INTEGER(I4B), INTENT(IN)                                :: jy

!  Internal variables

REAL(DP), PARAMETER                                     :: a1 = -0.004539  
REAL(DP), PARAMETER                                     :: a2 = -0.0001638  
REAL(DP), PARAMETER                                     :: a3 = 0.00002551
REAL(DP), PARAMETER                                     :: c1 = -9.9559    
REAL(DP), PARAMETER                                     :: c2 = 7.0845     
REAL(DP), PARAMETER                                     :: c3 = 3.9093
REAL(DP), PARAMETER                                     :: capa = -3.033405 
REAL(DP), PARAMETER                                     :: capb = 10.128163 
REAL(DP), PARAMETER                                     :: capc = -8.750567 
REAL(DP), PARAMETER                                     :: capd = 2.663107
REAL(DP)                                                :: tc
REAL(DP)                                                :: term1
REAL(DP)                                                :: term2
REAL(DP)                                                :: term3
REAL(DP)                                                :: cmolal
REAL(DP)                                                :: xsum

INTEGER(I4B)                                            :: ls

!  Data from Helgeson and Kirkham (1974), for 500 bars

tc = t(jx,jy,jz)

IF (DensityModule == 'temperature') THEN

! Data from Helgeson and Kirkham (1974), for 500 bars
!!  ro(jx,jy,jz) = 0.102116058D+04 - 0.494459726D-01*tc -  &
!!    0.486833333D-02*tc*tc + 0.143363552D-04*tc*tc*tc - 0.219537985D-07*(tc)**04
!!  ro(jx,jy,jz) = 997.0

!!McCutcheon, S.C., Martin, J.L, Barnwell, T.O. Jr. 1993. Water 
!!Quality in Maidment, D.R. (Editor). Handbood of Hydrology, 
!!McGraw-Hill, New York, NY (p. 11.3 )

  ro(jx,jy,jz) = 1000.0d0*(1.0d0 - (tc + 288.9414d0) / (508929.2d0*(tc + 68.12963d0))*(tc-3.9863d0)**2.0d0)

ELSE IF (DensityModule == 'sodium_chloride') THEN

  cmolal = SQRT( s(MeanSalt(1),jx,jy,jz)*s(MeanSalt(2),jx,jy,jz)  )
! Bethke's sodium chloride fit
  term1 = c1 * EXP(a1*cmolal)
  term2 = c2 * EXP(a2*tc)
  term3 = c3 * EXP(a3)          !  Assumes 1 bar pressure for now
!!  term3 = c3 * EXP(a3*pbar)   !  Original, giving dependence on pressure in bars
  xsum = term1 + term2 + term3
  ro(jx,jy,jz) = 1000.0*( capa + xsum * (capb + xsum * (capc + xsum * capd)) )

ELSE IF (DensityModule == 'sodium_nitrate') THEN                             !  Use linear regression based on molality from CRC
  cmolal = SQRT( s(MeanSalt(1),jx,jy,jz)*s(MeanSalt(2),jx,jy,jz)  )
  ro(jx,jy,jz) = 1000.0*(1.008974 + 0.04195*cmolal)
ELSE IF (DensityModule == 'potassium_nitrate') THEN                          !  Use linear regression based on molality from CRC
  cmolal = SQRT( s(MeanSalt(1),jx,jy,jz)*s(MeanSalt(2),jx,jy,jz)  )
  ro(jx,jy,jz) = 1000.0*(1.002226618 + 0.052782279*cmolal)
ELSE IF (DensityModule == 'calcium_nitrate') THEN                             !  Use linear regression based on molality from CRC
  cmolal = SQRT( 0.5*s(MeanSalt(1),jx,jy,jz)*s(MeanSalt(2),jx,jy,jz)  )
  ro(jx,jy,jz) = 1000.0*(1.008974 + 0.04195*cmolal)
ELSE IF (DensityModule == 'martin') THEN
  ro(jx,jy,jz) = 1000.0 + 25.0*s(MeanSalt(1),jx,jy,jz)*1000.0
ELSE
  CALL stringlen(DensityModule,ls)
  WRITE(*,*)
  WRITE(*,*) ' Density module not recognized: ', DensityModule(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

RETURN
END SUBROUTINE density
!  ***********************************************************
