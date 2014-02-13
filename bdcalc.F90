!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:53:00
 
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

SUBROUTINE bdcalc(ncomp,nspec,nbnd)
USE crunchtype
USE params
USE concentration
USE RunTime

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nbnd

! Internal arrays and variables

REAL(DP)                                      :: sum
REAL(DP)                                      :: MeanSaltConcentration

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: ksp

IF (DensityModule /= 'temperature') THEN
! Calculate the correction for the mass fraction of water:  kg_solution/kg_water 
  MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*sbnd(MeanSalt(1),nbnd) +   &
         wtaq(MeanSalt(2))*sbnd(MeanSalt(2),nbnd)) 
  xgrambnd(nbnd) = 1.0/(1.0 + MeanSaltConcentration)
ELSE
 xgrambnd(nbnd) = 1.0
END IF

DO i = 1,ncomp
  sum=0.0
  DO ksp = 1,nspec
    sum = sum + muaq(ksp,i)*spb(ksp+ncomp,nbnd)
  END DO
  s_local(i) = xgrambnd(nbnd)*(sum + spb(i,nbnd))
!!  s_local(i) = sum + spb(i,nbnd)
END DO

RETURN
END SUBROUTINE bdcalc
!**************************************************************
