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

SUBROUTINE read_constantgasflow(nout,nx,ny,nz,constant_gasflow,  &
    qxgasinit,qygasinit,qzgasinit)
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
LOGICAL(LGT), INTENT(OUT)                                   :: constant_gasflow
REAL(DP), INTENT(OUT)                                       :: qxgasinit
REAL(DP), INTENT(OUT)                                       :: qygasinit
REAL(DP), INTENT(OUT)                                       :: qzgasinit

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz
constant_gasflow = .false.
qxgasinit = 0.0
qygasinit = 0.0
qzgasinit = 0.0

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
  
!  Check to see if initial substring is "constant_gasflow"
  
  IF (ssch == 'constant_gasflow' .OR. ssch == 'constant_gas_flow') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      qxgasinit = 0.0
      qygasinit = 0.0
      qzgasinit = 0.0
      constant_gasflow = .true.
      RETURN
    END IF
    
!  First, looking for X velocity
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qxgasinit = DNUM(ssch)
        constant_gasflow = .true.
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "constant_gasflow" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "constant_gasflow" '
      WRITE(*,*) ' Value(s) must be given if constant_gasflow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
!  Now, looking for Y velocity if NY not equal to 1
    
    IF (ny == 1) THEN
      qygasinit = 0.0
      qzgasinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qygasinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following X gas velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Y gas velocity given'
      WRITE(*,*) ' Y velocity must be given if constant_gasflow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
! Now, looking for Z velocity if NZ not equal to 1
    
    IF (nz == 1) THEN
      qzgasinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qzgasinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following Y gas velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Z gas velocity given'
      WRITE(*,*) ' Z velocity must be given if constant_gasflow specified'
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
END SUBROUTINE read_constantgasflow
