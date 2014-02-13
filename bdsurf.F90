!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:53:09
 
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

SUBROUTINE bdsurf(ncomp,nsurf,nsurf_sec,nbnd,scorr)
USE crunchtype
USE concentration
USE params

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: nbnd
REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr

! Internal arrays and variables

REAL(DP)                                      :: sum

INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: ns

DO is = 1,nsurf
  sum=0.0D0
  DO ns = 1,nsurf_sec
    sum = sum + musurf(ns,is+ncomp)*spsurfb(ns+nsurf,nbnd)
  END DO
  scorr(is) = sum + spsurfb(is,nbnd)
END DO

RETURN
END SUBROUTINE bdsurf
!**************************************************************
