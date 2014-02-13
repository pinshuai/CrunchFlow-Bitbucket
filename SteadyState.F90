SUBROUTINE SteadyState(ncomp,nx,ny,nz,delt,steady)
USE crunchtype
USE concentration
USE runtime

IMPLICIT NONE

!   External variables and arrays

INTEGER(I4B), INTENT(IN)                                         :: ncomp
INTEGER(I4B), INTENT(IN)                                         :: nx
INTEGER(I4B), INTENT(IN)                                         :: ny
INTEGER(I4B), INTENT(IN)                                         :: nz
REAL(DP), INTENT(IN)                                             :: delt
LOGICAL(LGT), INTENT(OUT)                                        :: steady

!  Internal variables and arrays

INTEGER(I4B)                                                     :: jx
INTEGER(I4B)                                                     :: jy
INTEGER(I4B)                                                     :: jz
INTEGER(I4B)                                                     :: i
INTEGER(I4B)                                                     :: imax
INTEGER(I4B)                                                     :: jxmax
INTEGER(I4B)                                                     :: jymax
INTEGER(I4B)                                                     :: jzmax
INTEGER(I4B)                                                     :: jxwidth
INTEGER(I4B)                                                     :: jywidth
INTEGER(I4B)                                                     :: jzwidth
INTEGER(I4B)                                                     :: width
INTEGER(I4B)                                                     :: ls

REAL(DP)                                                         :: tderivative
REAL(DP)                                                         :: maxchange

maxchange = 0.0
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO i = 1,ncomp
         tderivative = (s(i,jx,jy,jz) - sn(i,jx,jy,jz))/(delt*sn(i,jx,jy,jz))
         IF (ABS(tderivative) > maxchange) THEN
           maxchange = ABS(tderivative)
           imax = i
           jxmax = jx
           jymax = jy
           jzmax = jz
         END IF
      END DO
    END DO
  END DO
END DO

CALL stringlen(ulab(imax),ls)

WRITE(*,*) 
WRITE(*,FMT=100) maxchange
WRITE(*,*) 'Component                                    :  ', ulab(imax)(1:ls)
WRITE(*,FMT=102) jxmax,jymax,jzmax
WRITE(*,*)

100 FORMAT(' Maximum rate of change (normalized mol/kg/yr): ',1pe12.4)
102 FORMAT(' Grid location                                : ',I3,':',I3,':',I3)

IF (maxchange < steadytol) THEN
  steady = .true.
ELSE
  steady = .false.
END IF

RETURN
END SUBROUTINE SteadyState
