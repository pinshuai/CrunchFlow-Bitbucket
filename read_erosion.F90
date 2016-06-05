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


SUBROUTINE read_erosion(nout,erodex,erodey)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), INTENT(OUT)                                       :: erodex
REAL(DP), INTENT(OUT)                                       :: erodey

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

REWIND nout

10  READ(nout,'(a)',END=1000) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ls == 0) THEN
  GO TO 10
END IF

lzs=ls
CALL convan(ssch,lzs,res)

!  Check to see if initial substring is "erode_x"

IF (ssch == 'erode_x') THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls == 0) THEN
    GO TO 10
  END IF
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    erodex = DNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' Cant interpret string following "erode_x" '
    WRITE(*,*) ' A numerical value should follow'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE IF (ssch == 'erode_y') THEN
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls == 0) THEN
    GO TO 10
  END IF
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    erodey = DNUM(ssch)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' Cant interpret string following "erode_y" '
    WRITE(*,*) ' A numerical value should follow'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  GO TO 10
END IF

GO TO 10

1000 RETURN
END SUBROUTINE read_erosion
