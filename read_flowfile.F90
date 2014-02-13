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

SUBROUTINE read_flowfile(nout,nx,ny,nz,constant_flow,  &
    qxinit,qyinit,qzinit,velocityfile,lfile,FileFormatType)
USE crunchtype
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
LOGICAL(LGT), INTENT(OUT)                                   :: constant_flow
REAL(DP), INTENT(OUT)                                       :: qxinit
REAL(DP), INTENT(OUT)                                       :: qyinit
REAL(DP), INTENT(OUT)                                       :: qzinit
CHARACTER (LEN=mls), INTENT(OUT)                            :: velocityfile
INTEGER(I4B), INTENT(OUT)                                   :: lfile
CHARACTER (LEN=mls), INTENT(OUT)                            :: FileFormatType

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat

nxyz = nx*ny*nz
velocityfile = ' '

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
  IF (ssch == 'read_velocity' .OR. ssch == 'read_velocityfile') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      qxinit = 0.0
      qyinit = 0.0
      qzinit = 0.0
      constant_flow = .true.
      RETURN
    END IF
    constant_flow = .false.
    
!  Looking for velocity file name (expecting prefix only)
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL stringtype(ssch,lzs,res)
      velocityfile = ssch
      lfile = ls

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
      WRITE(*,*) ' No file name following "read_VelocityFile" in FLOW section '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_flowfile
