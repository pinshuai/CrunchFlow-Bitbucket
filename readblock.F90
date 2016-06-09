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


SUBROUTINE readblock(nin,nout,section,found,ncount)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nin
INTEGER(I4B), INTENT(IN)                                    :: nout
CHARACTER (LEN=mls), INTENT(IN)                             :: section
LOGICAL(LGT), INTENT(OUT)                                   :: found
INTEGER(I4B), INTENT(OUT)                                   :: ncount

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2

INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nlen2
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs


dummy2 = ' '
found = .false.

REWIND nin
REWIND nout
ncount = 0

100 CONTINUE
READ(nin,'(a)',END=400) dummy1
nlen1 = LEN(dummy1)
nlen2 = LEN(section)
CALL majuscules(dummy1,nlen1)
IF (dummy1 == section) THEN
  found = .true.
  200   CONTINUE
  READ(nin,'(a)') dummy1
  BACKSPACE(nin)
  READ(nin,'(a)') dummy2
  id = 1
  iff = mls
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
  END IF
  IF (dummy2(1:1) /= '!' .AND. dummy1 /= 'end' .AND.  &
        dummy1 /= 'End' .AND. dummy1 /= 'END' .AND. dummy1 /= ' ') THEN
    WRITE(nout,'(a)') dummy2
    IF (res == 'a') THEN
      ncount = ncount + 1
    END IF
  ELSE IF (dummy1 == 'end' .OR. dummy1 == 'END'  &
        .OR. dummy1 == 'End') THEN
    WRITE(nout,*)
    WRITE(nout,*)
  END IF
  IF (dummy1 == 'end' .OR. dummy1 == 'END' .OR. dummy1 == 'End') THEN
    GO TO 300
  END IF
  GO TO 200
ELSE
  GO TO 100
END IF
300 CONTINUE

REWIND nin
REWIND nout

400 RETURN
END SUBROUTINE readblock
