!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:58:04
 
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

SUBROUTINE gases(ncomp,ngas,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE temperature
USE runtime, ONLY: Duan

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: ngas
INTEGER(I4B), INTENT(IN)                                   :: jx 
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz
!  Internal variables

REAL(DP)                                                   :: tempk
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: pg
REAL(DP)                                                   :: ln_fco2
REAL(DP)                                                   :: vrInOut

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: kk

ln_fco2 = 0.0d0

IF (Duan) THEN
  pg = GasPressureTotal(jx,jy,jz)
END IF

tempk = t(jx,jy,jz) + 273.15

!!denmol = LOG(1.e05/(8.314*tempk))   ! P/RT = n/V, with pressure converted from bars to Pascals
denmol = DLOG( (1.0E05) /(8.314d0*tempk) )   ! P/RT = n/V, with pressure converted from bars to Pascals

!!  NOTE:  The "denmol" should convert to mol/m*3 (n/V)

DO kk = 1,ngas

  IF (ikh2o /= 0) THEN

    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + mugas(kk,i)*(gam(i,jx,jy,jz))
      ELSE
         sum = sum + mugas(kk,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
      END IF
    END DO

  ELSE

    sum = 0.0
    DO i = 1,ncomp
      sum = sum + mugas(kk,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
    END DO

  END IF

  IF (Duan) THEN
    ln_fco2 = 0.0d0  ! fugacity coefficient for CO2(g)
    if (namg(kk) == 'CO2(g)') then
      vrInOut = vrSave(jx,jy,jz)
      call fugacity_co2(pg,tempk,ln_fco2,vrInOut)
      vrSave(jx,jy,jz) = vrInOut
    end if
  END IF

!! Basically, first two terms on RHS give you the mole fraction, then multipled by n/V to gives mol/m^3

  spgas(kk,jx,jy,jz) = keqgas(kk,jx,jy,jz) + sum + denmol - ln_fco2
  spgas10(kk,jx,jy,jz) = DEXP(spgas(kk,jx,jy,jz))  ! mol/m**3

END DO

RETURN
END SUBROUTINE gases
!  **************************************************************
