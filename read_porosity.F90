!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:07:17
 
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

SUBROUTINE read_porosity(nout,nco,portemp,porfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nco
REAL(DP), INTENT(OUT)                                       :: portemp
LOGICAL(LGT), INTENT(IN OUT)                                :: porfound

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: lcond

CHARACTER (LEN=mls)                                         :: dumstring


dumstring = condlabel(nco)
CALL stringlen(dumstring,lcond)

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'set_porosity' .OR. ssch == 'fix_porosity') THEN
    GO TO 200
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' "Set_porosity" should be followed by a numeric value'
    WRITE(*,*) '   In condition: ', dumstring(1:lcond)
    WRITE(*,*) '       ABORTING RUN  '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  porfound = .true.
  portemp = DNUM(ssch)
  RETURN
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No numerical value for porosity given '
  WRITE(*,*) '   In condition: ', dumstring(1:lcond)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

300 RETURN

END SUBROUTINE read_porosity
