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

SUBROUTINE readRecharge(nout,RechargeCondition)
USE crunchtype
USE params
USE concentration
USE transport
USE flow
USE strings
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: RechargeCondition


!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: intbnd_tmp

RechargeCondition = 0

REWIND nout

10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'recharge') THEN
      
!  Now, look for geochemical condition following infiltration
      
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
      WRITE(*,*) ' Geochemical condition for recharge in MODFLOW block not recognized'
      WRITE(*,*) ' Label: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
      50 intbnd_tmp = nco
    ELSE         !  Blank string
      WRITE(*,*)
      WRITE(*,*) ' No geochemical condition for recharge in MODFLOW block provided'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
      
    RechargeCondition = intbnd_tmp

  ELSE
    GO TO 10
  END IF
ELSE
  GO TO 10
END IF

500 RETURN
END SUBROUTINE readRecharge
