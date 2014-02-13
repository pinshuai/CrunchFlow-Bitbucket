SUBROUTINE readtypescan1(num)
! ***********************************************************************
!  read constant head, well, drain, or river results
! ***********************************************************************

USE crunchtype

IMPLICIT NONE

INTEGER(I4B), INTENT(OUT)                                :: num

CHARACTER (LEN=16)                                       :: text

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: i
INTEGER(I4B)                                             :: j
INTEGER(I4B)                                             :: k

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: idum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: jdum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: kdum
REAL(DP), DIMENSION(:), ALLOCATABLE                     :: yyy

BACKSPACE 1
READ(1) kper,kstp,ncol,nrow,nlay,text,num

ALLOCATE(idum(num))
ALLOCATE(jdum(num))
ALLOCATE(kdum(num))
ALLOCATE(yyy(num))

DO i =1,num
  READ(1) kdum(i),jdum(i),idum(i),yyy(i)
END DO

DEALLOCATE(idum)
DEALLOCATE(jdum)
DEALLOCATE(kdum)
DEALLOCATE(yyy)

RETURN
END SUBROUTINE readtypescan1

! ************* End SUBROUTINE readtypescan1 ********************************
