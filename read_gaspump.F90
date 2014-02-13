!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:06:55
 
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

SUBROUTINE read_gaspump(nout,nx,ny,nz,nchem,ngaspump)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: ngaspump

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: intbnd_tmp
INTEGER(I4B)                                                :: jxxtemp
INTEGER(I4B)                                                :: jyytemp
INTEGER(I4B)                                                :: jzztemp

REAL(DP)                                                    :: qtemp

nxyz = nx*ny*nz

REWIND nout

ngaspump = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'gaspump') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qtemp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "gaspump"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
!  Now, look for geochemical condition following pumping rate (only used if rate is positive)
      
      id = ids + ls
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
        WRITE(*,*) ' Geochemical condition for pumping well not found'
        WRITE(*,*) ' Label = ',ssch
        WRITE(*,*)
        READ(*,*)
        STOP
        50         ngaspump = ngaspump+ 1
        intbnd_tmp = nco
      ELSE         !  Blank string
        WRITE(*,*)
        WRITE(*,*) ' No geochemical condition for gas pumping well provided'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
! Now look for pumping well
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          jxxtemp = JNUM(ssch)
        ELSE                !  An ascii string--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' A grid location should follow gas pumping well specification'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE                  ! Zero length trailing string
        WRITE(*,*)
        WRITE(*,*) ' No grid location given for gas pumping zone'
        WRITE(*,*) ' Gas pumping zone ',ngaspump
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
      IF (jxxtemp > nx) THEN
        WRITE(*,*)
        WRITE(*,*) ' You have specified a gas pumping zone at JX > NX'
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        READ(*,*)
        STOP
      END IF
      IF (jxxtemp < 1) THEN
        WRITE(*,*)
        WRITE(*,*) ' You have specified a gas pumping zone at JX < 1'
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        READ(*,*)
        STOP
      END IF
      
      WRITE(*,*)
      WRITE(*,*) ' Gas pumping zone number ',ngaspump
      WRITE(*,*) ' Jxx location = ', jxxtemp
      WRITE(*,*)
      
!!      IF (ny > 1) THEN
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            jyytemp = JNUM(ssch)
          ELSE                !  An ascii string--so bag it.
            WRITE(*,*)
            WRITE(*,*) ' No Y location for gas pumping zone'
            WRITE(*,*) ' Gas pumping zone ',ngaspump
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE                  ! Zero length trailing string
          WRITE(*,*)
          WRITE(*,*) ' No Y location for gas pumping zone'
          WRITE(*,*) ' Gas pumping zone ',ngaspump
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
        IF (jyytemp > ny) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JY > NY'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        IF (jyytemp < 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JY < 1'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        
        WRITE(*,*)
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        WRITE(*,*) ' Jyy location = ', jyytemp
        WRITE(*,*)
        
!!      ELSE
!!        jyytemp = 1
!!      END IF


        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            jzztemp = JNUM(ssch)
          ELSE                !  An ascii string--so bag it.
            WRITE(*,*)
            WRITE(*,*) ' No Z location for gas pumping zone'
            WRITE(*,*) ' Gas pumping zone ',ngaspump
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE                  ! Zero length trailing string
          WRITE(*,*)
          WRITE(*,*) ' No Z location for gas pumping zone'
          WRITE(*,*) ' Gas pumping zone ',ngaspump
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
        IF (jzztemp > nz) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JZ > NZ'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        IF (jzztemp < 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JZ < 1'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        
        WRITE(*,*)
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        WRITE(*,*) ' Jzz location = ', jzztemp
        WRITE(*,*)

      gaspump(jxxtemp,jyytemp,jzztemp) = qtemp
      intbndgas(jxxtemp,jyytemp,jzztemp) = intbnd_tmp
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No gas pumping rate given'
      WRITE(*,*) ' Gas pumping zone ignored'
      WRITE(*,*)
      GO TO 10
    END IF
  ELSE
    GO TO 10
  END IF
ELSE
  GO TO 10
END IF

GO TO 10

500 RETURN
END SUBROUTINE read_gaspump
