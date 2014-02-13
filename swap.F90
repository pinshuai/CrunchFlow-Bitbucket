SUBROUTINE swap(a,b)
USE crunchtype
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(INOUT)                 :: a,b
REAL(DP), DIMENSION(SIZE(a))                          :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap
