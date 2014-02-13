!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:08
 
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

SUBROUTINE read_string(nout,l_string,parchar,parfind,dumstring,section)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(OUT)                                       :: l_string

CHARACTER (LEN=mls), INTENT(IN)                                 :: section
CHARACTER (LEN=mls), INTENT(IN)                                 :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                             :: parfind
CHARACTER (LEN=mls), INTENT(OUT)                                :: dumstring

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lchar

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
!      if(ls.ne.0) then
lzs=ls
CALL convan(ssch,lzs,res)
IF (ssch == parchar) THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
!      endif
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
l_string = ls
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a string'
    WRITE(*,*) ' In section ',section
    WRITE(*,*) ' Following parameter ',parchar
    WRITE(*,*) ' Aborting run'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  dumstring = ssch
  RETURN
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No trailing string found after'
  WRITE(*,*) parchar
  WRITE(*,*) ' In section ',section
  WRITE(*,*)
  parfind = ' '
END IF

RETURN
END SUBROUTINE read_string
