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
    
subroutine CrunchPETScInitializeDiffusion(nx,ny,nz,xvecD,bvecD,amatD,userD,ierr)
USE crunchtype
USE solver, ONLY:  XvecCrunchD, BvecCrunchD
USE transport, ONLY:  aDD,bDD,cDD,dDD,eDD,fDD,gDD,hDD,iDD

IMPLICIT NONE

#include "petsc/finclude/petsc.h"


!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

INTEGER(I4B)                                                          :: nxyz
INTEGER(I4B)                                                          :: linefil

! ******************** PETSC declarations ********************************
PetscFortranAddr     userD(*)
Mat                  amatD
Vec                  bvecD,xvecD
!!SLES                 sles
PC                   pc
KSP                  ksp
!!Scalar               zeroPetsc
! ************************end PETSc declarations of PETSc variables ******

nxyz = nx*ny*nz

IF (ALLOCATED(XvecCrunchD)) THEN
  DEALLOCATE(XvecCrunchD)
END IF
ALLOCATE(XvecCrunchD(0:nxyz-1))
IF (ALLOCATED(BvecCrunchD)) THEN
  DEALLOCATE(BvecCrunchD)
END IF
ALLOCATE(BvecCrunchD(0:nxyz-1))

IF (ny == 1 .AND. nz == 1) THEN       !  1D problem (assumes an X coordinate direction
  linefil = 3
ELSE
  IF (ny > 1 .AND. nz > 1) THEN      !  3D problem
    linefil = 7
  ELSE
    linefil = 5
  END IF
END IF

call MatCreateSeqAIJ(PETSC_COMM_SELF,nxyz,nxyz,linefil,PETSC_NULL_INTEGER,amatD,ierr)
!!call MatSetOption(amatD,MAT_COLUMN_ORIENTED,ierr) 
call MatSetOption(amatD,MAT_ROW_ORIENTED,PETSC_FALSE,ierr) 
call VecCreateSeqWithArray(PETSC_COMM_SELF,1,nxyz,BvecCrunchD,bvecD,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,1,nxyz,XvecCrunchD,xvecD,ierr)
call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
call KSPSetFromOptions(ksp,ierr)
call KSPGetPC(ksp,pc,ierr)
!!call SLESGetKSP(sles,ksp,ierr)

userD(1) = amatD
userD(2) = bvecD
userD(3) = xvecD
!!userD(4) = sles
userD(5) = pc
userD(6) = ksp

  IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN            ! 1D case
    IF (ALLOCATED(aDD)) THEN
      DEALLOCATE(aDD)
    END IF
    ALLOCATE(aDD(nx,ny,nz))
    IF (ALLOCATED(bDD)) THEN
      DEALLOCATE(bDD)
    END IF
    ALLOCATE(bDD(nx,ny,nz))
    IF (ALLOCATED(cDD)) THEN
      DEALLOCATE(cDD)
    END IF
    ALLOCATE(cDD(nx,ny,nz))
  ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN         ! 3D case
    IF (ALLOCATED(aDD)) THEN
      DEALLOCATE(aDD)
    END IF
    ALLOCATE(aDD(nx,ny,nz))
    IF (ALLOCATED(bDD)) THEN
      DEALLOCATE(bDD)
    END IF
    ALLOCATE(bDD(nx,ny,nz))
    IF (ALLOCATED(cDD)) THEN
      DEALLOCATE(cDD)
    END IF
    ALLOCATE(cDD(nx,ny,nz))
    IF (ALLOCATED(dDD)) THEN
      DEALLOCATE(dDD)
    END IF
    ALLOCATE(dDD(nx,ny,nz))
    IF (ALLOCATED(eDD)) THEN
      DEALLOCATE(eDD)
    END IF
    ALLOCATE(eDD(nx,ny,nz))
    IF (ALLOCATED(fDD)) THEN
      DEALLOCATE(fDD)
    END IF
    ALLOCATE(fDD(nx,ny,nz))
    IF (ALLOCATED(gDD)) THEN
      DEALLOCATE(gDD)
    END IF
    ALLOCATE(gDD(nx,ny,nz))
    IF (ALLOCATED(hDD)) THEN
      DEALLOCATE(hDD)
    END IF
    ALLOCATE(hDD(nx,ny,nz))
    IF (ALLOCATED(iDD)) THEN
      DEALLOCATE(iDD)
    END IF
    ALLOCATE(iDD(nx,ny,nz))
  ELSE
    IF (ny > 1 .AND. nz == 1) THEN
      IF (ALLOCATED(aDD)) THEN
        DEALLOCATE(aDD)
      END IF
      ALLOCATE(aDD(nx,ny,nz))
      IF (ALLOCATED(bDD)) THEN
        DEALLOCATE(bDD)
      END IF
      ALLOCATE(bDD(nx,ny,nz))
      IF (ALLOCATED(cDD)) THEN
        DEALLOCATE(cDD)
      END IF
      ALLOCATE(cDD(nx,ny,nz))
      IF (ALLOCATED(dDD)) THEN
        DEALLOCATE(dDD)
      END IF
      ALLOCATE(dDD(nx,ny,nz))
      IF (ALLOCATED(eDD)) THEN
        DEALLOCATE(eDD)
      END IF
      ALLOCATE(eDD(nx,ny,nz))
      IF (ALLOCATED(fDD)) THEN
        DEALLOCATE(fDD)
      END IF
      ALLOCATE(fDD(nx,ny,nz))
    ELSE
      IF (ALLOCATED(aDD)) THEN
        DEALLOCATE(aDD)
      END IF
      ALLOCATE(aDD(nx,ny,nz))
      IF (ALLOCATED(bDD)) THEN
        DEALLOCATE(bDD)
      END IF
      ALLOCATE(bDD(nx,ny,nz))
      IF (ALLOCATED(cDD)) THEN
        DEALLOCATE(cDD)
      END IF
      ALLOCATE(cDD(nx,ny,nz))
      IF (ALLOCATED(gDD)) THEN
        DEALLOCATE(gDD)
      END IF
      ALLOCATE(gDD(nx,ny,nz))
      IF (ALLOCATED(hDD)) THEN
        DEALLOCATE(hDD)
      END IF
      ALLOCATE(hDD(nx,ny,nz))
      IF (ALLOCATED(iDD)) THEN
        DEALLOCATE(iDD)
      END IF
      ALLOCATE(iDD(nx,ny,nz))
    END IF
  END IF

RETURN
END subroutine CrunchPETScInitializeDiffusion
