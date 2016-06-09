!! CrunchTope 
!! Copyright (c) 2016, Carl Steefel
!! Copyright (c) 2016, The Regents of the University of California, 
!! through Lawrence Berkeley National Laboratory (subject to 
!! receipt of any required approvals from the U.S. Dept. of Energy).  
!! All rights reserved.

!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met: 

!! (1) Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.

!! (2) Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.

!! (3) Neither the name of the University of California, Lawrence
!! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of    
!! its contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.

!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE   
    
subroutine CrunchPETScTolerances(user,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
USE crunchtype
USE solver, ONLY:  level, SolverMethod, PCMethod
IMPLICIT NONE

#include "petsc/finclude/petsc.h"

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

