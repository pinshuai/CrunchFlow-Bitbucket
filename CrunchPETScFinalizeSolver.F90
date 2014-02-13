subroutine CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
USE crunchtype

IMPLICIT NONE

#include "finclude/petsc.h"

!  External variables and arrays

INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays


! ******************** PETSC declarations ********************************
PetscFortranAddr     userC(*)
Vec                  xvec,bvec
Mat                  amatpetsc
KSP                  ksp
! ************************end PETSc declarations of PETSc variables ******


ksp = userC(6)

call VecDestroy(xvec,ierr)
call VecDestroy(bvec,ierr)
call MatDestroy(amatpetsc,ierr)
call KSPDestroy(ksp,ierr)

RETURN
END subroutine CrunchPETScFinalizeSolver
