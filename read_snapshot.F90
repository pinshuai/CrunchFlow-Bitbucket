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


SUBROUTINE read_snapshot(nout,lchar,parchar,parchar2,parfind,realmult,lenarray,section)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: lchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar2
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
INTEGER(I4B), INTENT(OUT)                                   :: lenarray

REAl(DP), DIMENSION(:), INTENT(IN OUT)                      :: realmult

CHARACTER (LEN=mls), INTENT(IN)                             :: section

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: npar

REAL(DP)                                                    :: realjunk

LOGICAL(LGT)                                                :: continuation

continuation = .FALSE.

REWIND nout

npar = 0

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
ELSE
  GOTO 100
END IF
IF (ssch == parchar .OR. ssch == parchar2) THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 CONTINUE

IF (continuation) THEN
  id = 1
  continuation = .false.
ELSE
  id = ids + ls
END IF


CALL sschaine(zone,id,iff,ssch,ids,ls)

IF (ssch == '&') THEN
  READ(nout,'(a)') zone
  continuation = .TRUE.
  GOTO 200
END IF

IF(ls /= 0) THEN
  npar = npar + 1
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a numeric value'
    WRITE(*,*) ' In section ',section
    WRITE(*,*) ' Following parameter ',parchar
    WRITE(*,*) ' Aborting run'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  realmult(npar) = DNUM(ssch)
  lenarray = npar
  GO TO 200
ELSE
  npar = npar - 1
  realjunk = -500.0
END IF

RETURN
END SUBROUTINE read_snapshot
