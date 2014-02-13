!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:45
 
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

SUBROUTINE totsurf_init(ncomp,nsurf,nsurf_sec)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec

!  Internal variables and arrays

INTEGER(I4B)                                                :: is
INTEGER(I4B)                                                :: ns

REAL(DP)                                                    :: sum

DO is = 1,nsurf
  sum=0.0
  DO ns = 1,nsurf_sec
    sum = sum + musurf(ns,is+ncomp)*spsurftmp10(ns+nsurf)
  END DO
  ssurftmp(is) = sum + spsurftmp10(is)
  continue
END DO

RETURN
END SUBROUTINE totsurf_init
!  **************************************************************
