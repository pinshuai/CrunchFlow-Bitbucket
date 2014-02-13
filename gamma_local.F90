!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:52
 
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

SUBROUTINE gamma_local(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec
INTEGER(I4B), INTENT(IN)                                   :: jx 
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz

!  Internal variables

REAL(DP)                                                   :: ctotal
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: sion_tmp
REAL(DP)                                                   :: tempc

INTEGER(I4B)                                               :: ik

ctotal = 0.0
DO ik = 1,ncomp+nspec
  ctotal = ctotal + sp10(ik,jx,jy,jz)
END DO

tempc = t(jx,jy,jz)

IF (ntemp == 1) THEN
  ah = adh(1)
  bh = bdh(1)
  bdt = bdot(1)
ELSE
  ah = adhcoeff(1) + adhcoeff(2)*tempc  &
      + adhcoeff(3)*tempc*tempc + adhcoeff(4)*tempc*tempc*tempc  &
      + adhcoeff(5)*tempc*tempc*tempc*tempc
  bh = bdhcoeff(1) + bdhcoeff(2)*tempc  &
      + bdhcoeff(3)*tempc*tempc + bdhcoeff(4)*tempc*tempc*tempc  &
      + bdhcoeff(5)*tempc*tempc*tempc*tempc
  bdt = bdtcoeff(1) + bdtcoeff(2)*tempc  &
      + bdtcoeff(3)*tempc*tempc + bdtcoeff(4)*tempc*tempc*tempc  &
      + bdtcoeff(5)*tempc*tempc*tempc*tempc
END IF

sum = 0.0D0

DO ik = 1,ncomp+nspec
  sum = sum + sp10(ik,jx,jy,jz)*chg(ik)*chg(ik)
END DO
sion_tmp = 0.50D0*sum

DO ik = 1,ncomp+nspec
  IF (chg(ik) == 0.0) THEN
!!    gam_local(ik) = clg*0.1*sion_tmp
    gam_local(ik) = 0.0
    IF (ulab(ik) == 'H2O') THEN
      gam_local(ik) = LOG(1.0/ctotal)
    END IF
  ELSE
    gam_local(ik) = -clg*( (ah*chg(ik)*chg(ik)*SQRT(sion_tmp))/  &
        (1.0+ acmp(ik)*bh*SQRT(sion_tmp)) + bdt*sion_tmp )
  END IF
END DO

RETURN
END SUBROUTINE gamma_local
!*****************************************************************
