!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:06:59
 
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

SUBROUTINE read_series(nout,nseries,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: nseries
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: i
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: workchar1

IF (ALLOCATED(TimeSeriesFile)) THEN
  DEALLOCATE(TimeSeriesFile)
END IF
ALLOCATE(TimeSeriesFile(1))

nxyz = nx*ny*nz

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
  
!  Check to see if initial substring is "time_series" (could be more than 1)
  
  IF (ssch == 'time_series' .OR. ssch == 'time_series_at_node') THEN
    nseries = nseries + 1
  ELSE
    GO TO 10
  END IF
ELSE         ! No string found
  GO TO 10
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  i = size(TimeSeriesFile,1)
  ALLOCATE(workchar1(i))
  workchar1 = TimeSeriesFile
  DEALLOCATE(TimeSeriesFile)
  ALLOCATE(TimeSeriesFile(nseries))
  IF(nseries /= 0) TimeSeriesFile(1:nseries-1) = workchar1(1:nseries-1)
  DEALLOCATE(workchar1)
  TimeSeriesFile(nseries) = ssch
ELSE
  IF (nseries == 1) THEN
    TimeSeriesFile(1) = 'TimeSeries.out'
  ELSE
    WRITE(*,*) 
    WRITE(*,*) ' File name for time series required when more than one is used'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (nxyz == 1) THEN
  jxseries(nseries) = 1
  jyseries(nseries) = 1
  jzseries(nseries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    jxseries(nseries) = JNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location (X) should follow "time_series" label'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (ny == 1) THEN
  jyseries(nseries) = 1
  jzseries(nseries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    jyseries(nseries) = JNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location (Y) should follow "time_series" label'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Y location for timeseries must be specified'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (nz == 1) THEN
  jzseries(nseries) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    jzseries(nseries) = JNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location (Z) should follow "time_series" label'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Z location for timeseries must be specified'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF 

500 GO TO 10

1000  RETURN
END SUBROUTINE read_series
