!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:32
 
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

SUBROUTINE read_condition(nin,nout,found,ncount,nchem,endoffile)
USE crunchtype
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nin
INTEGER(I4B), INTENT(IN)                                    :: nout

LOGICAL(LGT), INTENT(OUT)                                   :: found
INTEGER(I4B), INTENT(OUT)                                   :: ncount
INTEGER(I4B), INTENT(IN)                                    :: nchem
LOGICAL(LGT), INTENT(IN OUT)                                :: endoffile

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1

dummy2 = ' '
found = .false.

REWIND nout
ncount = 0

100 CONTINUE
READ(nin,'(a)',END=400) dummy1
nlen1 = LEN(dummy1)
CALL majuscules(dummy1,nlen1)
!  First, parse the string (dummy1) to see if first word
!  is "condition"
id = 1
iff = mls
CALL sschaine(dummy1,id,iff,ssch,ids,ls)
IF (ssch == 'condition') THEN
  found = .true.
! Check for label following "condition"
  id = ids + ls
!  WRITE(*,*) ' Checking for geochemical condition label'
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    condlabel(nchem) = ssch
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No label for condition in input file'
    WRITE(*,5050) nchem
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
! Now, get the "condition" title (use whatever part of the main
!   string is left over)
  
  id = ids + ls
  condtitle(nchem) = dummy1(id:mls)
  
! Now, read the contents of the condition block
  
  200   CONTINUE       ! Return here while reading condition block
  READ(nin,'(a)') dummy1
  BACKSPACE(nin)
  READ(nin,'(a)') dummy2
  
  IF (dummy2(1:1) /= '!' .AND. dummy1 /= 'end' .AND.  &
        dummy1 /= 'End' .AND. dummy1 /= 'END') THEN
    WRITE(nout,'(a)') dummy2
    ncount = ncount + 1
  ELSE IF (dummy1 == 'end' .OR. dummy1 == 'END'  &
        .OR. dummy1 == 'End') THEN
    WRITE(nout,*)
    WRITE(nout,*)
  END IF
  IF (dummy1 == 'end' .OR. dummy1 == 'END' .OR. dummy1 == 'End') THEN
    GO TO 300    ! Exit loop
  END IF
  GO TO 200
ELSE
  GO TO 100
END IF

300 CONTINUE    ! Exit here

REWIND nout
RETURN

400 endoffile = .true.

5050 FORMAT(1X,'Condition number ',i2,' in input file')

RETURN
END SUBROUTINE read_condition
