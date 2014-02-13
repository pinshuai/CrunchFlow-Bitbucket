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

SUBROUTINE readModFlowStress(nout,perlen,nstp,tsmult,ndimdummy)
USE crunchtype
USE CrunchFunctions
USE strings
USE modflowModule

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: perlen
INTEGER(I4B), DIMENSION(:), INTENT(OUT)                     :: nstp
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: tsmult
INTEGER(I4B), INTENT(IN)                                    :: ndimdummy

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: i

REWIND nout

nstress = 0

10  READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ssch == 'modflow_stress') THEN
  nstress = nstress + 1
  IF (nstress > ndimdummy) THEN
    WRITE(*,*)
    WRITE(*,*) ' Too many stress periods'
    WRITE(*,*) ' Redimension "ndimdummy" in routine "readModFlowStress" '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Read the length of this stress period

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      perlen(nstress) = DNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading length of stress period in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading length of stress period in MODFLOW block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Now, read the number of time steps in the stress period

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      nstp(nstress) = JNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading number of time steps in stress period in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading number of time steps in stress period'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Now, read the multiplier

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      tsmult(nstress) = DNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading time step multiplier in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading time step multiplier in MODFLOW block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

END IF

GO TO 10

500  RETURN

END SUBROUTINE readModFlowStress
