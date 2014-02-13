!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:50
 
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

SUBROUTINE jacobian_plus(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver
USE transport

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz

!  Internal variables

REAL(DP)                                             :: sum1
REAL(DP)                                             :: sum2
REAL(DP)                                             :: sum3
REAL(DP)                                             :: tmp
REAL(DP)                                             :: spec_conc
REAL(DP)                                             :: mutemp
REAL(DP)                                             :: dcorrectmp
REAL(DP)                                             :: chgtmp

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp

fjac(:,:,jx,jy,jz) = 0.0
fjac_d(:,:,jx,jy,jz) = 0.0
fjac_chg(:,:,jx,jy,jz) = 0.0

DO ksp = 1,nspec
  spec_conc = sp10(ksp+ncomp,jx,jy,jz)
  dcorrectmp = d_correct(ksp+ncomp)
  chgtmp = chg(ksp+ncomp)
  DO i = 1,ncomp
    IF (muaq(ksp,i) /= 0.0) THEN
      mutemp = muaq(ksp,i)
      DO i2 = 1,ncomp
        fjac(i2,i,jx,jy,jz) = fjac(i2,i,jx,jy,jz) +   &
              mutemp*muaq(ksp,i2)*spec_conc
        fjac_d(i2,i,jx,jy,jz) = fjac_d(i2,i,jx,jy,jz) +   &
              dcorrectmp*mutemp*muaq(ksp,i2)*spec_conc
        fjac_chg(i2,i,jx,jy,jz) = fjac_chg(i2,i,jx,jy,jz) +   &
              chgtmp*dcorrectmp*mutemp*muaq(ksp,i2)*spec_conc
      END DO
    END IF
  END DO
END DO

DO i = 1,ncomp
  tmp = sp10(i,jx,jy,jz)
  fjac(i,i,jx,jy,jz) = fjac(i,i,jx,jy,jz) + tmp
  fjac_d(i,i,jx,jy,jz) = fjac_d(i,i,jx,jy,jz) + tmp*d_correct(i)
  fjac_chg(i,i,jx,jy,jz) = fjac_chg(i,i,jx,jy,jz) + tmp*d_correct(i)*chg(i)
END DO

RETURN
END SUBROUTINE jacobian_plus
!***********************************************************
