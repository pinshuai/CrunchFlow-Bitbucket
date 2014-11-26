subroutine CrunchPETScTolerances(user,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
USE crunchtype
USE solver, ONLY:  level, SolverMethod, PCMethod
IMPLICIT NONE

#include "finclude/petsc.h"

!  External variables and arrays
REAL(DP), INTENT(IN)                                                  :: rtolksp
REAL(DP), INTENT(IN)                                                  :: atolksp
REAL(DP), INTENT(IN)                                                  :: dtolksp

INTEGER(I4B), INTENT(IN)                                              :: maxitsksp
INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

! ******************** PETSC declarations ********************************
PetscFortranAddr     user(*)
!!SLES                 sles
PC                   pc
KSP                  ksp
! ************************end PETSc declarations of PETSc variables ******

!!sles = user(4)
pc = user(5)
ksp = user(6)
 
! Tolerances for linear solver set here

!!!call KSPSetTolerances(ksp,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
call KSPSetTolerances(ksp,rtolksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,maxitsksp,ierr)

! Choose linear solver

IF (SolverMethod == 'gmres') THEN
  CALL KSPSetType(ksp,KSPGMRES,ierr)
ELSE IF (SolverMethod == 'bcgs') THEN
  CALL KSPSetType(ksp,KSPBCGS,ierr)
ELSE IF (SolverMethod == 'direct') THEN
  CALL KSPSetType(ksp,KSPPREONLY,ierr)
ELSE 
  CALL KSPSetType(ksp,KSPBICG,ierr)
END IF

! Choose preconditioning method

IF (PCMethod == 'jacobi') THEN
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCJACOBI,ierr)
  CALL PCFactorSetLevels(pc,level,ierr)
ELSE
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCILU,ierr)
  CALL PCFactorSetLevels(pc,level,ierr)
END IF

IF (SolverMethod == 'direct') THEN
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCCHOLESKY,ierr)
END IF

RETURN
END subroutine CrunchPETScTolerances

