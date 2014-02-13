SUBROUTINE AllocateGasesOS3D(nx,ny,nz,ncomp)

USE crunchtype
USE params
USE runtime
USE concentration
USE transport
USE solver
USE ReadFlow

IMPLICIT NONE

!!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: nx
INTEGER(I4B), INTENT(IN)                               :: ny
INTEGER(I4B), INTENT(IN)                               :: nz
INTEGER(I4B), INTENT(IN)                               :: ncomp

IF (ALLOCATED(sgas)) THEN
  DEALLOCATE(sgas)
END IF
ALLOCATE(sgas(ncomp,nx,ny,nz))

IF (ALLOCATED(sgasn)) THEN
  DEALLOCATE(sgasn)
END IF
ALLOCATE(sgasn(ncomp,nx,ny,nz))

IF (ALLOCATED(fgas_local)) THEN
  DEALLOCATE(fgas_local)
END IF
ALLOCATE(fgas_local(ncomp,ncomp))

sgas = 0.0
sgasn = 0.0
fgas_local = 0.0

RETURN
END SUBROUTINE AllocateGasesOS3D