!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:04:53
 
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

SUBROUTINE readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
              conditionNum,modflowstring,lenstring)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: nparams
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jxTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jyTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jzTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: conditionNum
CHARACTER (LEN=mls), INTENT(IN)                             :: modflowstring
INTEGER(I4B), INTENT(IN)                                    :: lenstring


!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
REWIND nout

nparams = 0
conditionNum = 0

10  READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ssch == modflowstring) THEN
  nparams = nparams + 1

  !  Read the location (JX:JY:JZ) following the "well" keyword

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      jxTemp(nparams) = JNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
      WRITE(*,*) '     A JX location should follow keyword ',modflowstring(1:lenstring)
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
    WRITE(*,*) '     A JX location should follow keyword ',modflowstring(1:lenstring)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      jyTemp(nparams) = JNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
      WRITE(*,*) '     A JY location should follow keyword ',modflowstring(1:lenstring)
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
      WRITE(*,*) '     A JY location should follow keyword ',modflowstring(1:lenstring)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      jzTemp(nparams) = JNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
      WRITE(*,*) '     A JZ location should follow keyword ',modflowstring(1:lenstring)
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
    WRITE(*,*) '     A JZ location should follow keyword ',modflowstring(1:lenstring)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
 
!  Now read the geochemical condition

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)

    IF (ssch == 'notneeded') THEN
      conditionNum(nparams) = 1
    ELSE

!     Check to see that geochemical condition is recognized

      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          conditionNum(nparams) = nco
          GO TO 550
        END IF
      END DO

      WRITE(*,*)
      WRITE(*,*) ' Geochemical condition label for parameter not recognized: ',modflowstring(1:lenstring)
      WRITE(*,*) '   Label: ',ssch(1:ls)
      WRITE(*,501) jxTemp(nparams),jyTemp(nparams),jzTemp(nparams)
      WRITE(*,*)
      READ(*,*)
      STOP

    END IF

  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading parameters in MODFLOW block'
    WRITE(*,*) '     No label for the geochemical condition found following MODFLOW parameter: ', &
                       modflowstring(1:lenstring)
    WRITE(*,501) jxTemp(nparams),jyTemp(nparams),jzTemp(nparams)

    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  550 CONTINUE

END IF

GO TO 10

501 FORMAT(1x,' MODFLOW parameter location: ',i4,1x,i4,1x,i4)

500  RETURN
END SUBROUTINE readModFlowParameters
