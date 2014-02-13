!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:53:03
 
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

SUBROUTINE bdexchange(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nbnd,sexb)
USE crunchtype
USE concentration
USE params

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN OUT)                  :: nbnd
REAL(DP), DIMENSION(:), INTENT(OUT)           :: sexb

!  Internal variables and arrays

REAL(DP)                                      :: sum

INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: nex
INTEGER(I4B)                                  :: ns

DO i = 1,ncomp
  sum=0.0
  DO nex = 1,nexch_sec
    sum = sum + muexc(nex,i)*spexb(nex+nexchange,nbnd)
  END DO
  DO ns = 1,nsurf_sec
    sum = sum + musurf(ns,i)*spsurfb(ns+nsurf,nbnd)
  END DO
  sexb(i) = sum
END DO

RETURN
END SUBROUTINE bdexchange
!**************************************************************
