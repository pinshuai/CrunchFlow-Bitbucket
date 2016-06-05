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

    
SUBROUTINE PauseFor(PauseLength)
USE crunchtype

IMPLICIT NONE

!  External variables and arrays

REAL(DP), INTENT(IN)                                       :: PauseLength

!  Internal variables and arrays

CHARACTER (LEN=12)                                         :: dumm1
CHARACTER (LEN=12)                                         :: dumm2
CHARACTER (LEN=12)                                         :: dumm3
INTEGER(I4B), DIMENSION(8)                                 :: curr_time

INTEGER(I4B)                                               :: str_mon
INTEGER(I4B)                                               :: str_day
INTEGER(I4B)                                               :: str_hr
INTEGER(I4B)                                               :: str_min
INTEGER(I4B)                                               :: str_sec
INTEGER(I4B)                                               :: end_mon
INTEGER(I4B)                                               :: end_day
INTEGER(I4B)                                               :: end_hr
INTEGER(I4B)                                               :: end_min
INTEGER(I4B)                                               :: end_sec

REAL(DP)                                                   :: simu_t

curr_time = 0
CALL date_and_time(dumm1,dumm2,dumm3,curr_time)
str_mon = curr_time(2)
str_day = curr_time(3)
str_hr  = curr_time(5)
str_min = curr_time(6)
str_sec = curr_time(7)

100 CALL date_and_time(dumm1,dumm2,dumm3,curr_time)
end_mon = curr_time(2)
end_day = curr_time(3)
end_hr  = curr_time(5)
end_min = curr_time(6)
end_sec = curr_time(7)
           
simu_t=(end_day*24*3600 + end_hr*3600 + end_min*60 + end_sec) -   &
                 (str_day*24*3600 + str_hr*3600 + str_min*60 + str_sec)

IF (simu_t < PauseLength) THEN
  GOTO 100
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE PauseFor
