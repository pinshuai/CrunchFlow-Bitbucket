!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:08:28
 
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

SUBROUTINE squeeze(substring,ls)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

CHARACTER (LEN=mls), INTENT(IN OUT)                           :: substring
INTEGER(I4B), INTENT(OUT)                                     :: ls

CHARACTER (LEN=mls)                                           :: dummy

INTEGER(I4B)                                                  :: i

ls=0
DO i=1,mls
  IF(substring(i:i) == ' ') THEN
    ls=ls +1
  ELSE
    ls = ls+1
    dummy = substring(ls:mls)
    substring = dummy
!          write(*,*) substring
!          pause
    GO TO 100
  END IF
END DO

100 RETURN
END SUBROUTINE squeeze

