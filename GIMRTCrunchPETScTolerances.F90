subroutine GIMRTCrunchPETScTolerances(user,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
USE crunchtype
USE solver, ONLY:  GIMRTlevel, GIMRT_SolverMethod, GIMRT_PCMethod
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

call KSPSetTolerances(ksp,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

! Choose linear solver

IF (GIMRT_SolverMethod == 'bcgs') THEN
  CALL KSPSetType(ksp,KSPBCGS,ierr)
ELSE 
  CALL KSPSetType(ksp,KSPGMRES,ierr)
END IF

! Choose preconditioning method

IF (GIMRT_PCMethod == 'ilu') THEN
  CALL PCSetType(pc,PCILU,ierr)
  CALL PCFactorSetLevels(pc,GIMRTlevel,ierr)
ELSE
  CALL PCSetType(pc,PCBJACOBI,ierr)
  CALL PCFactorSetLevels(pc,GIMRTlevel,ierr)
END IF

RETURN
END subroutine GIMRTCrunchPETScTolerances

