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
