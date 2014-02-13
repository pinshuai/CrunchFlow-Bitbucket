SUBROUTINE readtype3(nx,ny,nz,text,xxx)
! ***********************************************************************
!  read thickness and X,Y,Z flows
! ***********************************************************************

USE crunchtype

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

CHARACTER (LEN=16), INTENT(OUT)                          :: text

REAL(DP), DIMENSION(nx,ny,nz), INTENT(OUT)                :: xxx

INTEGER(I4B) :: griddims,temp

REAL(SIP), DIMENSION(nx,ny,nz)                            :: yyy

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: i
INTEGER(I4B)                                             :: j
INTEGER(I4B)                                             :: k


!  Add 1 to gridDims to account for satThick

BACKSPACE 1
READ(1,END = 1000) kper,kstp,ncol,nrow,nlay,text
READ(1)yyy

IF (text == 'THKSAT' .AND. yyy(1,1,1) == -111.0) THEN
  CONTINUE
ELSE
  xxx=yyy
END IF

RETURN
1000  CONTINUE
END SUBROUTINE readtype3

! ************* End SUBROUTINE readtype3 *********************************
