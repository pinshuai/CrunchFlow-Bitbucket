subroutine CrunchPETScInitializePressure(nx,ny,nz,userP,ierr,xvecP,bvecP,amatP)
USE crunchtype
USE flow, ONLY:  XvecCrunchP, BvecCrunchP

IMPLICIT NONE

#include "finclude/petsc.h"

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

INTEGER(I4B)                                                          :: nxyz
INTEGER(I4B)                                                          :: linefil

! ******************** PETSC declarations ********************************
PetscFortranAddr     userP(*)
Mat                  amatP
Vec                  bvecP,xvecP
!!SLES                 sles
PC                   pc
KSP                  ksp
!!Scalar               zeroPetsc
! ************************end PETSc declarations of PETSc variables ******

nxyz = nx*ny*nz

IF (ALLOCATED(XvecCrunchP)) THEN
  DEALLOCATE(XvecCrunchP)
END IF
ALLOCATE(XvecCrunchP(0:nxyz-1))
IF (ALLOCATED(BvecCrunchP)) THEN
  DEALLOCATE(BvecCrunchP)
END IF
ALLOCATE(BvecCrunchP(0:nxyz-1))

IF (ny == 1 .AND. nz == 1) THEN       !  1D problem (assumes an X coordinate direction
  linefil = 3
ELSE
  IF (ny > 1 .AND. nz > 1) THEN      !  3D problem
    linefil = 7
  ELSE
    linefil = 5
  END IF
END IF

call MatCreateSeqAIJ(PETSC_COMM_SELF,nxyz,nxyz,linefil,PETSC_NULL_INTEGER,amatP,ierr)
!!call MatSetOption(amatP,MAT_COLUMN_ORIENTED,ierr)
call MatSetOption(amatP,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,nxyz,BvecCrunchP,bvecP,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,nxyz,XvecCrunchP,xvecP,ierr)
call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
call KSPSetFromOptions(ksp,ierr)
call KSPGetPC(ksp,pc,ierr)
!!call SLESGetKSP(sles,ksp,ierr)

userP(1) = amatP
userP(2) = bvecP
userP(3) = xvecP
!!userP(4) = sles
userP(5) = pc
userP(6) = ksp

RETURN
END subroutine CrunchPETScInitializePressure
