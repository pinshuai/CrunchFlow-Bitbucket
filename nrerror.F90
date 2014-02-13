SUBROUTINE nrerror(string)
USE crunchtype
IMPLICIT NONE
CHARACTER (LEN=*), intent(in)                 :: string
WRITE (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'

END SUBROUTINE nrerror
