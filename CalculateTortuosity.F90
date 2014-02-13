SUBROUTINE CalculateTortuosity(nx,ny,nz)
USE crunchtype
USE params
USE runtime
USE crunch_interface
USE medium
USE transport

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                                           :: nx
INTEGER(I4B), INTENT(IN)                                           :: ny
INTEGER(I4B), INTENT(IN)                                           :: nz

!!  Internal arrays and variables

INTEGER(I4B)                                                       :: ls
INTEGER(I4B)                                                       :: jx
INTEGER(I4B)                                                       :: jy
INTEGER(I4B)                                                       :: jz

REAL(DP)                                                           :: PercolationThreshold
REAL(DP)                                                           :: Coefficient
REAL(DP)                                                           :: LimitingTortuosity

!!        *****************************

IF (TortuosityOption == 'CostaRica' .OR. TortuosityOption == 'costarica' .OR. TortuosityOption == 'COSTARICA') THEN

  PercolationThreshold = 0.16d0
  Coefficient = 0.2d0
  LimitingTortuosity = 2.0D-05

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        IF (por(jx,jy,jz) >= PercolationThreshold) THEN
          tortuosity(jx,jy,jz) = LimitingTortuosity + Coefficient*(por(jx,jy,jz) - PercolationThreshold)**2.0
        ELSE
          tortuosity(jx,jy,jz) = LimitingTortuosity
        END IF
      END DO
    END DO
  END DO
ELSE
  CALL stringlen(TortuosityOption,ls)
  WRITE(*,*)
  WRITE(*,*) ' Tortuosity option not recognized: ',TortuosityOption(1:ls)
  WRITE(*,*)
  STOP
END IF

RETURN

END SUBROUTINE CalculateTortuosity