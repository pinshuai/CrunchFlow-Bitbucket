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

SUBROUTINE readFileName(nout,l_string,parchar,parfind,dumstring,section,FileFormatType)
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
CHARACTER (LEN=mls), INTENT(OUT)                                :: FileFormatType

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lchar
INTEGER(I4B)                                                    :: lensection
INTEGER(I4B)                                                    :: lenparchar
INTEGER(I4B)                                                    :: lenformat

CALL stringlen(section,lensection)
CALL stringlen(parchar,lenparchar)

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
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a file name'
    WRITE(*,*) '   In section: ',section(1:lensection)
    WRITE(*,*) '   Following parameter: ',parchar(1:lenparchar)
    WRITE(*,*) ' Aborting run'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  dumstring = ssch

!!  Now, check for a file format
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  lenformat = ls
  CALL majuscules(ssch,lenformat)
  IF (ls /= 0) THEN
    IF (ssch == 'singlecolumn') THEN
      FileFormatType = 'SingleColumn'
    ELSE IF (ssch == 'continuousread') THEN
      FileFormatType = 'ContinuousRead'
    ELSE IF (ssch == 'unformatted') THEN
      FileFormatType = 'Unformatted' 
    ELSE IF (ssch == 'distanceplusvariable' .OR. ssch == 'fullform' .OR. ssch == 'full') THEN
      FileFormatType = 'FullForm'         
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' File format not recognized: ', ssch(1:lenformat)
      WRITE(*,*)
      READ(*,*)
      STOP
    ENDIf
  ELSE    !! No file format provided, so assume default
    FileFormatType = 'SingleColumn'
  ENDIF
  RETURN
ELSE
  WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a file name'
    WRITE(*,*) '   In section: ',section(1:lensection)
    WRITE(*,*) '   Following parameter: ',parchar(1:lenparchar)
  WRITE(*,*)
  parfind = ' '
END IF

RETURN
END SUBROUTINE readFileName
