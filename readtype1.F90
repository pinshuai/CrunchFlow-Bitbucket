SUBROUTINE readtype1(ndim,iii,jjj,kkk,xxx,num)
! ***********************************************************************
!  read constant head, well, drain, or river results
! ***********************************************************************

USE crunchtype

IMPLICIT NONE

! External arrays and variables
INTEGER(I4B), INTENT(IN)                               :: ndim
INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: iii
INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: jjj
INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: kkk

REAL(DP), DIMENSION(:), INTENT(OUT)                    :: xxx
INTEGER(I4B), INTENT(OUT)                              :: num

!  Internal arrays and variables

CHARACTER (LEN=16)                                       :: text

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: i
INTEGER(I4B)                                             :: j
INTEGER(I4B)                                             :: k

REAL(SIP), DIMENSION(:), ALLOCATABLE                     :: yyy

BACKSPACE 1
READ(1) kper,kstp,ncol,nrow,nlay,text,num

ALLOCATE(yyy(num))

IF (num > ndim) THEN
  WRITE(*,*)
  WRITE(*,*) ' Arrays in "readtype1" routine not dimensioned large enough'
  WRITE(*,*) ' Current dimension: ', ndim
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

DO i =1,num
  READ(1) kkk(i),jjj(i),iii(i),yyy(i)
END DO

DO i = 1,num
  xxx(i) = yyy(i)
END DO

DEALLOCATE(yyy)

RETURN
END SUBROUTINE readtype1

! ************* End SUBROUTINE readtype1 ********************************
