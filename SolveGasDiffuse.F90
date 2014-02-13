SUBROUTINE SolveGasDiffuse(nx,ny,nz,nn,icomp,delt,user,amatD)
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
REAL(DP)                                                              :: portemp
REAL(DP)                                                              :: satgas

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

IF (icomp > 1) GOTO 500               !  No need to calculate matrix, since matrix has been calculated for component 1

call MatZeroEntries(amatD,ierr)

IF (nx > 1 .AND. ny ==1 .AND. nz == 1) THEN           ! 1D problem assuming jx is coordinate

  jy = 1
  jz = 1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1  
    satgas = 1.0 - satliq(jx,jy,jz)
    portemp = por(jx,jy,jz)    
    IF (activecell(jx,jy,jz) == 0) THEN
      DiagonalTerm = 1.0
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
    ELSE
      AccumulationTerm = dxy(jx,jy,jz)*portemp*satgas/delt
      DiagonalTerm = bg(jx,jy,jz) + AccumulationTerm
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,ag(jx,jy,jz),INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,cg(jx,jy,jz),INSERT_VALUES,ierr)
    END IF
  END DO

  jx = 1
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  satgas = 1.0 - satliq(jx,jy,jz)
  portemp = por(jx,jy,jz)   
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
  ELSE
    AccumulationTerm = dxy(jx,jy,jz)*portemp*satgas/delt
    DiagonalTerm = bg(jx,jy,jz) + AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,cg(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

  jx = nx
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  satgas = 1.0 - satliq(jx,jy,jz)
  portemp = por(jx,jy,jz)
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
  ELSE
    AccumulationTerm = dxy(jx,jy,jz)*portemp*satgas/delt
    DiagonalTerm = bg(jx,jy,jz) + AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,ag(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

ELSE                                                !  2D problem

    jz = 1
    DO jy = 1,ny
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-1,ag(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,cg(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jx = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+1,0.0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+1,cg(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jx = nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-1,0.0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-1,ag(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    jz = 1
    DO jx = 1,nx
      DO jy = 2,ny-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-nx,0.0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,0.0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-nx,fg(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,dg(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jy = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+nx,0.0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+nx,dg(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jy = ny
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-nx,0.0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-nx,fg(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          satgas = 1.0 - satliq(jx,jy,jz)
          portemp = por(jx,jy,jz) 
          j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
          IF (activecell(jx,jy,jz) == 0) THEN
            DiagonalTerm = 1.0
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          ELSE
            AccumulationTerm = dxy(jx,jy,jz)*portemp*satgas/delt
            DiagonalTerm = bg(jx,jy,jz) + eg(jx,jy,jz) + AccumulationTerm
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          END IF
        END DO
      END DO
    END DO

END IF

500 CONTINUE

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      satgas = 1.0 - satliq(jx,jy,jz)
      portemp = por(jx,jy,jz) 
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
      IF (activecell(jx,jy,jz) == 0) THEN
        RightHandSide = sgasn(icomp,jx,jy,jz)
      ELSE
        RightHandSide = dxy(jx,jy,jz)*portemp*satgas*sgasn(icomp,jx,jy,jz)/delt
      END IF
      BvecCrunchD(j) = RightHandSide
    END DO
  END DO
END DO

CALL MatAssemblyBegin(amatD,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatD,MAT_FINAL_ASSEMBLY,ierr)

RETURN
END SUBROUTINE SolveGasDiffuse
