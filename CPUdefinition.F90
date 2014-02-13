subroutine def_time(CPU_unit)

USE crunchtype

IMPLICIT NONE

REAL(DP),INTENT(OUT) :: CPU_unit
REAL(DP)                :: T1,T2
REAL(DP),DIMENSION(:,:), ALLOCATABLE :: A, B, C

ALLOCATE(A(1000,1000),B(1000,1000),C(1000,1000))

call random_number(A)
call random_number(B)

call CPU_time(T1)
C=MATMUL(A,B)
call cpu_time(T2)

PRINT*,'T1,T2',T1,T2
CPU_unit=T2-T1
PRINT*,'CPU_unit',CPU_unit

DEALLOCATE (A)
DEALLOCATE (B)
DEALLOCATE (C)

return
end subroutine  def_time