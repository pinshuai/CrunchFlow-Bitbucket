!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:09:23
 
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

SUBROUTINE totexchange(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nexchange
INTEGER(I4B), INTENT(IN)                                    :: nexch_sec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: nex
INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: ix

REAL(DP)                                                    :: sum

DO i = 1,ncomp
  sum=0.0
  DO nex = 1,nexch_sec
    sum = sum + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
  END DO
  DO ns = 1,nsurf_sec
    sum = sum + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
  END DO
  sch(i,jx,jy,jz) = sum
END DO

      do ix = 1,nexchange
        sum = 0.0
        do nex = 1,nexch_sec
          sum = sum + muexc(nex,ix+ncomp)*spex10(nex+nexchange,jx,jy,jz)
        end do
        totex(ix) = sum
      end do

RETURN
END SUBROUTINE totexchange
!  **************************************************************
