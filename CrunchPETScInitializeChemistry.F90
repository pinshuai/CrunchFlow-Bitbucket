subroutine CrunchPETScInitializeChemistry(nx,ny,nz,neqn,xvec,bvec,amatpetsc,userC,ierr)
USE crunchtype
USE solver, ONLY:  xn,fxx

IMPLICIT NONE

#include "finclude/petsc.h"

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN)                                              :: neqn

INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

INTEGER(I4B)                                                          :: linefil
INTEGER(I4B)                                                          :: np
INTEGER(I4B)       :: i
INTEGER(I4B), Dimension(nx)      :: iRowWidth

! ******************** PETSC declarations ********************************
PetscFortranAddr     userC(*)
Mat                  amatpetsc
Vec                  bvec,xvec
!!SLES                 sles
PC                   pc
KSP                  ksp
!!Scalar               zeroPetsc
! ************************end PETSc declarations of PETSc variables ******

np = neqn*nx*ny*nz

IF (ny == 1 .AND. nz == 1) THEN       !  1D problem (assumes an X coordinate direction
  linefil = 3
ELSE
  IF (ny > 1 .AND. nz > 1) THEN      !  3D problem
    linefil = 7
    WRITE(*,*)
    WRITE(*,*) ' 3D not yet allowed for GIMRT calculation'
    WRITE(*,*)
    call PetscFinalize(ierr)
    STOP
  ELSE
    linefil = 5
  END IF
END IF

!!!do i = 2,nx-1
!!!  iRowWidth(i) = 3
!!!end do
!!!iRowWidth(1) = 2
!!!iRowWidth(nx) = 2

call MatCreateSeqBAIJ(PETSC_COMM_SELF,neqn,np,np,linefil,PETSC_NULL_INTEGER,amatpetsc,ierr)
call MatSetOption(amatpetsc,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,neqn,np,fxx,bvec,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,neqn,np,xn,xvec,ierr)
call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
call KSPSetFromOptions(ksp,ierr)
call KSPGetPC(ksp,pc,ierr)


userC(1) = amatpetsc
userC(2) = bvec
userC(3) = xvec
!!userC(4) = sles
userC(5) = pc
userC(6) = ksp

RETURN
END subroutine CrunchPETScInitializeChemistry
