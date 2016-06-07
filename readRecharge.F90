!! CrunchTope 
!! Copyright (c) 2016, Carl Steefel
!! Copyright (c) 2016, The Regents of the University of California, 
!! through Lawrence Berkeley National Laboratory (subject to 
!! receipt of any required approvals from the U.S. Dept. of Energy).  
!! All rights reserved.

!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met: 

!! (1) Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.

!! (2) Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.

!! (3) Neither the name of the University of California, Lawrence
!! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of    
!! its contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.

!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE 


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
