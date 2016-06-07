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

    
SUBROUTINE PestScaleOutput(nout)
USE crunchtype
USE CrunchFunctions
USE params
USE strings
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: idum

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  IF (ssch == 'PestScaleOutput' .OR. ssch == 'pestscaleoutput') THEN
    GO TO 200
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  DO idum = 1,nplot
    IF (ssch == ulab(iplot(idum))) THEN
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL stringtype(ssch,lzs,res)
        IF (res /= 'n') THEN
          WRITE(*,*)
          WRITE(*,*) '  Species name in PestScaleOutput keyword should be followed by a numeric value'
          WRITE(*,*) '       ABORTING RUN  '
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        PestScale(idum)= DNUM(ssch)
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No numerical value for PestScaleOutput given'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END IF
  END DO
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No information given after PestScaleOutput keyword'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

GOTO 100

300 RETURN

END SUBROUTINE PestScaleOutput