SUBROUTINE readtype2(nx,ny,ll,qq)
! ***********************************************************************
!  read recharge or evapotranspiration results
! ***********************************************************************

USE crunchtype

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), DIMENSION(nx,ny), INTENT(INOUT)            :: ll


REAL(DP), DIMENSION(nx,ny), INTENT(INOUT)                :: qq
  
CHARACTER (LEN=16)                                       :: text

REAL(SIP), DIMENSION(nx,ny)                               :: rr

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: i
INTEGER(I4B)                                             :: j
INTEGER(I4B)                                             :: k
INTEGER(I2B), DIMENSION(nx,ny)                           :: ldum


BACKSPACE 1
READ(1) kper,kstp,ncol,nrow,nlay,text

READ(1) ((ldum(i,j),i=1,nx),j=1,ny)

READ(1) ((rr(i,j),i=1,nx),j=1,ny)

qq=rr
ll = ldum

RETURN
END SUBROUTINE readtype2

! ************* End SUBROUTINE readtype2 ********************************

