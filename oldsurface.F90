!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:01:45
 
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

SUBROUTINE oldsurface(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz

!  Internal variables and arrays

REAL(DP)                                      :: sum2
REAL(DP)                                      :: sum3

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: ksp
INTEGER(I4B)                                  :: nex
INTEGER(I4B)                                  :: ns

DO i = 1,ncomp
  sum2 = 0.0
  sum3 = 0.0
  DO nex = 1,nexch_sec
    sum2 = sum2 + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
  END DO
  DO ns = 1,nsurf_sec
    sum3 = sum3 + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
  END DO
  sexold(i,jx,jy,jz) = sum2
  ssurfold(i,jx,jy,jz) = sum3
END DO

RETURN
END SUBROUTINE oldsurface
!  **************************************************************
