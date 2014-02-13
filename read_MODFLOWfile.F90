!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:03:23
 
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

SUBROUTINE read_MODFLOWfile(nout,lfile,mxwell,mxrivr,mxdrn)
USE crunchtype
USE params
USE strings
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout

INTEGER(I4B), INTENT(OUT)                                   :: lfile
INTEGER(I4B), INTENT(OUT)                                   :: mxwell
INTEGER(I4B), INTENT(OUT)                                   :: mxrivr
INTEGER(I4B), INTENT(OUT)                                   :: mxdrn

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

CHARACTER (LEN=mls)                                         :: file

LOGICAL(LGT)                                                :: ext

MODFLOWfile = ' '
modflow = .FALSE.

REWIND nout

10  READ(nout,'(a)',END=1000) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res) 
  IF (ssch == 'read_modflow') THEN

!  Switch on MODFLOW option only if file found

    modflow = .TRUE.
    
!  Looking for MODFLOW file name (expecting prefix only)
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      MODFLOWfile = ssch
      lfile = ls
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name following "read_modflow" in MODFLOW section '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
    GO TO 10
  END IF
  
ELSE         ! No string found
  GO TO 10
END IF

IF (MODFLOWfile == ' ') THEN
  RETURN
ELSE
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".hff"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*) 
    WRITE(*,*) ' MODFLOW "*.hff" file not found: ', file(1:lfile+4)
    WRITE(*,*) 
    READ(*,*)
    STOP
  ELSE
    OPEN(1,file=file,form='unformatted',status='old')
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".bas"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*) 
    WRITE(*,*) ' MODFLOW "*.bas" file not found: ', file(1:lfile+4)
    WRITE(*,*)
    READ(*,*) 
    STOP
  ELSE
    OPEN(52,file=file,status='old')
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".wel"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxwell = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxwell
    CLOSE(53) 
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".drn"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxdrn = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxdrn
    CLOSE(53) 
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".riv"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxrivr = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxrivr
    CLOSE(53) 
  END IF

!!  file = ' '
!!  file(1:lfile)=MODFLOWfile(1:lfile)
!!  file(lfile+1:lfile+1+6)=".hffout"
!!  INQUIRE(FILE=file,EXIST=ext)
!!  IF (.NOT. ext) THEN
!!    WRITE(*,*) 
!!    WRITE(*,*) ' MODFLOW "*.hffout" file not found: ', file(1:lfile+7)
!!    WRITE(*,*) 
!!    STOP
!!  ELSE
!!    OPEN(3,file=file)
!!  END IF
END IF

1000 RETURN
END SUBROUTINE read_MODFLOWfile

