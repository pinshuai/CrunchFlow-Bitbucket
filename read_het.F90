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

SUBROUTINE read_het(nout,nchem,nhet,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: nhet
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
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b

nxyz = nx*ny*nz

REWIND nout

nhet = 0
10  READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check to see that heterogeneity label matches one of the labels
!  for geochemical conditions (condlabel)
  
  DO nco = 1,nchem
    IF (ssch == condlabel(nco)) THEN
      GO TO 50
    END IF
  END DO
  WRITE(*,*)
  WRITE(*,*) ' Label for heterogeneity not found in list of condition labels'
  WRITE(*,*) ' Label = ',ssch(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
  50   nhet = nhet + 1
  ndist(nhet) = nco
  IF (nxyz == 1) THEN      ! If NXYZ = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'fix') THEN
        jjfix(nhet) = 1
      ELSE
        jjfix(nhet) = 0
      END IF
    END IF
    jxxlo(nhet) = 1
    jxxhi(nhet) = 1
    jyylo(nhet) = 1
    jyyhi(nhet) = 1
    jzzlo(nhet) = 1
    jzzhi(nhet) = 1
    RETURN
  END IF
ELSE
  GO TO 10
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jxxlo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jxxhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jxxhi(nhet) = jxxlo(nhet)   !  Assume jxxhi=jxxlo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF


IF (ny == 1 .AND. nz == 1) THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch == 'fix') THEN
      jjfix(nhet) = 1
    ELSE
      jjfix(nhet) = 0
    END IF
  END IF
  jyylo(nhet) = 1
  jyyhi(nhet) = 1
  jzzlo(nhet) = 1
  jzzhi(nhet) = 1
  GO TO 10
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jyylo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jyyhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jyyhi(nhet) = jyylo(nhet)   !  Assume jyyhi=jyylo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (nz == 1) THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch == 'fix') THEN
      jjfix(nhet) = 1
    ELSE
      jjfix(nhet) = 0
    END IF
  END IF
  jzzlo(nhet) = 1
  jzzhi(nhet) = 1
  GO TO 10
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jzzlo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jzzhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jzzhi(nhet) = jzzlo(nhet)   !  Assume jzzhi=jzzlo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'fix') THEN
    jjfix(nhet) = 1
  ELSE
    jjfix(nhet) = 0
  END IF
END IF

GO TO 10

500  RETURN
END SUBROUTINE read_het
