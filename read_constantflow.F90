!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:35
 
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

SUBROUTINE read_constantflow(nout,nx,ny,nz,constant_flow,qxinit,qyinit,qzinit)
USE crunchtype
USE CrunchFunctions
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

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz
constant_flow = .false.
qxinit = 0.0
qyinit = 0.0
qzinit = 0.0

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check to see if initial substring is "constant_flow"
  
  IF (ssch == 'constant_flow' .OR. ssch == 'constant_liquidflow') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      qxinit = 0.0
      qyinit = 0.0
      qzinit = 0.0
      constant_flow = .true.
      RETURN
    END IF
    
!  First, looking for X velocity
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qxinit = DNUM(ssch)
        constant_flow = .true.
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "constant_flow" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "constant_flow" '
      WRITE(*,*) ' Value(s) must be given if constant_flow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
!  Now, looking for Y velocity if NY not equal to 1
    
    IF (ny == 1) THEN
      qyinit = 0.0
      qzinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qyinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following X velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Y velocity given'
      WRITE(*,*) ' Y velocity must be given if constant_flow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
! Now, looking for Z velocity if NZ not equal to 1
    
    IF (nz == 1) THEN
      qzinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qzinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following Y velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Z velocity given'
      WRITE(*,*) ' Z velocity must be given if constant_flow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    100     CONTINUE
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_constantflow
