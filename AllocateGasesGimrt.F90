SUBROUTINE AllocateGasesGimrt(nx,ny,nz,ncomp)

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

IF (ALLOCATED(ag)) THEN
  DEALLOCATE(ag)
END IF
ALLOCATE(ag(nx,ny,nz))

IF (ALLOCATED(bg)) THEN
  DEALLOCATE(bg)
END IF
ALLOCATE(bg(nx,ny,nz))

IF (ALLOCATED(cg)) THEN
  DEALLOCATE(cg)
END IF
ALLOCATE(cg(nx,ny,nz))

IF (ALLOCATED(dg)) THEN
  DEALLOCATE(dg)
END IF
ALLOCATE(dg(nx,ny,nz))

IF (ALLOCATED(eg)) THEN
  DEALLOCATE(eg)
END IF
ALLOCATE(eg(nx,ny,nz))

IF (ALLOCATED(fg)) THEN
  DEALLOCATE(fg)
END IF
ALLOCATE(fg(nx,ny,nz))

ag = 0.0
bg = 0.0
cg = 0.0
dg = 0.0
eg = 0.0
fg = 0.0

IF (ALLOCATED(sgas)) THEN
  DEALLOCATE(sgas)
END IF
ALLOCATE(sgas(ncomp,nx,ny,nz))

IF (ALLOCATED(sgasn)) THEN
  DEALLOCATE(sgasn)
END IF
ALLOCATE(sgasn(ncomp,nx,ny,nz))

IF (ALLOCATED(fgas)) THEN
  DEALLOCATE(fgas)
END IF
ALLOCATE(fgas(ncomp,ncomp,nx,ny,nz))

    sgas = 0.0
    sgasn = 0.0
    fgas = 0.0

RETURN
END SUBROUTINE AllocateGasesGimrt