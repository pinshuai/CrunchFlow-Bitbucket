!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:19
 
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

SUBROUTINE read_temperature(nout,FileFormatType)
USE crunchtype
USE params
USE temperature
USE strings

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                        :: nout
CHARACTER (LEN=mls), INTENT(OUT)                                :: FileFormatType

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lenformat

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'read_temperature' .OR. ssch == 'read_temperaturefile') THEN
    GO TO 200
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  tfile = ssch
  WRITE(*,*)
  WRITE(*,*) ' Temperature to be read from: ',ssch(1:ls)
  WRITE(*,*)
!!  *****************************************************************
!!    Now, check for a file format
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
!!  *****************************************************************
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Read temperature specified, but no file name provided'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

RETURN
END SUBROUTINE read_temperature
