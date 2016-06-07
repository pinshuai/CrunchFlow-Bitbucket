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

SUBROUTINE read_condition(nin,nout,found,ncount,nchem,endoffile)
USE crunchtype
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nin
INTEGER(I4B), INTENT(IN)                                    :: nout

LOGICAL(LGT), INTENT(OUT)                                   :: found
INTEGER(I4B), INTENT(OUT)                                   :: ncount
INTEGER(I4B), INTENT(IN)                                    :: nchem
LOGICAL(LGT), INTENT(IN OUT)                                :: endoffile

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dummy1
CHARACTER (LEN=mls)                                         :: dummy2
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1

dummy2 = ' '
found = .false.

REWIND nout
ncount = 0

100 CONTINUE
READ(nin,'(a)',END=400) dummy1
nlen1 = LEN(dummy1)
CALL majuscules(dummy1,nlen1)
!  First, parse the string (dummy1) to see if first word
!  is "condition"
id = 1
iff = mls
CALL sschaine(dummy1,id,iff,ssch,ids,ls)
IF (ssch == 'condition') THEN
  found = .true.
! Check for label following "condition"
  id = ids + ls
!  WRITE(*,*) ' Checking for geochemical condition label'
  CALL sschaine(dummy1,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    condlabel(nchem) = ssch
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No label for condition in input file'
    WRITE(*,5050) nchem
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
! Now, get the "condition" title (use whatever part of the main
!   string is left over)
  
  id = ids + ls
  condtitle(nchem) = dummy1(id:mls)
  
! Now, read the contents of the condition block
  
  200   CONTINUE       ! Return here while reading condition block
  READ(nin,'(a)') dummy1
  BACKSPACE(nin)
  READ(nin,'(a)') dummy2
  
  IF (dummy2(1:1) /= '!' .AND. dummy1 /= 'end' .AND.  &
        dummy1 /= 'End' .AND. dummy1 /= 'END') THEN
    WRITE(nout,'(a)') dummy2
    ncount = ncount + 1
  ELSE IF (dummy1 == 'end' .OR. dummy1 == 'END'  &
        .OR. dummy1 == 'End') THEN
    WRITE(nout,*)
    WRITE(nout,*)
  END IF
  IF (dummy1 == 'end' .OR. dummy1 == 'END' .OR. dummy1 == 'End') THEN
    GO TO 300    ! Exit loop
  END IF
  GO TO 200
ELSE
  GO TO 100
END IF

300 CONTINUE    ! Exit here

REWIND nout
RETURN

400 endoffile = .true.

5050 FORMAT(1X,'Condition number ',i2,' in input file')

RETURN
END SUBROUTINE read_condition
