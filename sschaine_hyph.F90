!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:08:42
 
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

SUBROUTINE sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
USE crunchtype
USE params
 
IMPLICIT NONE

!  External variables and arrays

CHARACTER (LEN=mls), INTENT(OUT)                           :: ssch_a
CHARACTER (LEN=mls), INTENT(OUT)                           :: ssch_b
CHARACTER (LEN=mls), INTENT(IN)                            :: zone

INTEGER(I4B), INTENT(IN)                                   :: id
INTEGER(I4B), INTENT(IN)                                   :: iff

INTEGER(I4B), INTENT(OUT)                                  :: ids
INTEGER(I4B), INTENT(OUT)                                  :: ls
INTEGER(I4B), INTENT(OUT)                                  :: ls_a
INTEGER(I4B), INTENT(OUT)                                  :: ls_b


!  Internal variables and arrays

CHARACTER (LEN=1)                                          :: delim
CHARACTER (LEN=1)                                          :: delim2

INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: i

ls = 0
ls_a=0
ls_b=0
ssch_a=' '
ssch_b=' '
DO i=id,iff
  IF(zone(i:i) /= ' ') THEN
    is=i
    ls = 1
    ls_a=1
    ls_b = 1
    ids=i
    delim='-'
    delim2=' '
    DO WHILE(zone(is:is) /= delim)
      IF (zone(is:is) == delim2) THEN
        ls=ls-1
        ls_a=ls_a-1
        ls_b = 0
        RETURN
      END IF
      ssch_a(ls_a:ls_a)=zone(is:is)
      ls_a=ls_a+1
      ls=ls+1
      is=is+1
    END DO
    ls = ls + 1
    ls_a = ls_a - 1
    IF (zone(is:is) == delim) THEN
      is=is+1
      delim = ' '
      DO WHILE(zone(is:is) /= delim)
        ssch_b(ls_b:ls_b)=zone(is:is)
        ls=ls+1
        ls_b=ls_b+1
        is=is+1
      END DO
      ls_b = ls_b - 1
      ls = ls - 1
    END IF
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE sschaine_hyph

