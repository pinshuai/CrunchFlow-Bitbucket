SUBROUTINE SolveDiffuseSpherical(nx,ny,nz,nn,icomp,delt,user,amatD)
USE crunchtype
USE params
USE concentration
USE solver
USE medium
USE transport
USE temperature

IMPLICIT NONE

!*****************************PETSc include statements ********************

#include "finclude/petsc.h"

!**************************** End PETSc include statements **************

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN)                                              :: nn
INTEGER(I4B), INTENT(IN)                                              :: icomp

REAL(DP), INTENT(IN)                                                  :: delt

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz
INTEGER(I4B)                                                          :: j
INTEGER(I4B)                                                          :: i
INTEGER(I4B)                                                          :: ierr
INTEGER(I4B)                                                          :: itsiterate
INTEGER(I4B)                                                          :: nxyz
     
REAL(DP)                                                              :: AccumulationTerm
REAL(DP)                                                              :: RightHandSide
REAL(DP)                                                              :: DiagonalTerm


! *******************begin PETSc declarations of f90 variables***********
INTEGER(I4B)             ::numprocs
INTEGER(I4B)             ::irank
INTEGER(I4B)             ::linefil
INTEGER(I4B), PARAMETER  ::maxitsksp=100
REAL(DP), PARAMETER      ::zero=0.0

!*********************end PETSc declarations ******************************

! ******************** PETSC declarations ********************************
PetscFortranAddr     user(6)
Mat                  amatD
! ************************end PETSc declarations of PETSc variables ******


IF (nn == 0) THEN
  RETURN
END IF

IF (icomp > 1) GOTO 500               !  No need to calculate matrix

call MatZeroEntries(amatD,ierr)


  jy = 1
  jz = 1
  DO jx = 2,nx-1
    j = jx-1   
    IF (activecell(jx,jy,jz) == 0) THEN
      DiagonalTerm = 1.0
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
    ELSE    
      AccumulationTerm = ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
      DiagonalTerm = b(jx,jy,jz) + dxy(jx,jy,jz)*AccumulationTerm
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,a(jx,jy,jz),INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,c(jx,jy,jz),INSERT_VALUES,ierr)
    END IF
  END DO

  jx = 1
  j = jx - 1
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
  ELSE 
    AccumulationTerm = ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
    DiagonalTerm = b(jx,jy,jz) + dxy(jx,jy,jz)*AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,c(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

  jx = nx
  j = jx - 1
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
  ELSE
    AccumulationTerm = ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
    DiagonalTerm = b(jx,jy,jz) + dxy(jx,jy,jz)*AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,a(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

500 CONTINUE

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
      IF (activecell(jx,jy,jz) == 0) THEN
        RightHandSide = sn(icomp,jx,jy,jz)
      ELSE
        RightHandSide = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)*sn(icomp,jx,jy,jz)/delt
      END IF
      BvecCrunchD(j) = RightHandSide
    END DO
  END DO
END DO

CALL MatAssemblyBegin(amatD,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatD,MAT_FINAL_ASSEMBLY,ierr)

RETURN
END SUBROUTINE SolveDiffuseSpherical
