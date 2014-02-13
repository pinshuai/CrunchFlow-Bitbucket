!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:45
 
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

SUBROUTINE readblock(nin,nout,section,found,ncount)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nin
INTEGER(I4B), INTENT(IN)                                    :: nout
CHARACTER (LEN=mls), INTENT(IN)                             :: section
LOGICAL(LGT), INTENT(OUT)                                   :: found
INTEGER(I4B), INTENT(OUT)                                   :: ncount

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2

INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nlen2
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs


dummy2 = ' '
found = .false.

REWIND nin
REWIND nout
ncount = 0

100 CONTINUE
READ(nin,'(a)',END=400) dummy1
nlen1 = LEN(dummy1)
nlen2 = LEN(section)
CALL majuscules(dummy1,nlen1)
IF (dummy1 == section) THEN
  found = .true.
  200   CONTINUE
  READ(nin,'(a)') dummy1
  BACKSPACE(nin)
  READ(nin,'(a)') dummy2
  id = 1
  iff = mls
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
  END IF
  IF (dummy2(1:1) /= '!' .AND. dummy1 /= 'end' .AND.  &
        dummy1 /= 'End' .AND. dummy1 /= 'END' .AND. dummy1 /= ' ') THEN
    WRITE(nout,'(a)') dummy2
    IF (res == 'a') THEN
      ncount = ncount + 1
    END IF
  ELSE IF (dummy1 == 'end' .OR. dummy1 == 'END'  &
        .OR. dummy1 == 'End') THEN
    WRITE(nout,*)
    WRITE(nout,*)
  END IF
  IF (dummy1 == 'end' .OR. dummy1 == 'END' .OR. dummy1 == 'End') THEN
    GO TO 300
  END IF
  GO TO 200
ELSE
  GO TO 100
END IF
300 CONTINUE

REWIND nin
REWIND nout

400 RETURN
END SUBROUTINE readblock
