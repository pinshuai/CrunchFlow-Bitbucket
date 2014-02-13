SUBROUTINE AllocateALL(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
USE crunchtype
USE params
USE runtime
USE crunch_interface
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE ReadFlow
USE modflowModule
USE NanoCrystal

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: npot
INTEGER(I4B), INTENT(IN)                           :: neqn
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz

!!  Internal arrays and variables

IF (ALLOCATED(sion)) THEN
  DEALLOCATE(sion)
END IF
ALLOCATE(sion(nx,ny,nz))
sion = 0.0
IF (ALLOCATED(fxmax)) THEN
  DEALLOCATE(fxmax)
END IF
ALLOCATE(fxmax(neqn))
fxmax = 0.0
IF (ALLOCATED(sumrd)) THEN
  DEALLOCATE(sumrd)
END IF
ALLOCATE(sumrd(ncomp+nexchange+nsurf))
sumrd = 0.0
IF (ALLOCATED(sumjackin)) THEN
  DEALLOCATE(sumjackin)
END IF
ALLOCATE(sumjackin(ncomp))
sumjackin = 0.0
IF (ALLOCATED(surf)) THEN
  DEALLOCATE(surf)
END IF
ALLOCATE(surf(nrct))
surf = 0.0
IF (ALLOCATED(jac_sat)) THEN
  DEALLOCATE(jac_sat)
END IF
ALLOCATE(jac_sat(ncomp+nexchange+nsurf))
jac_sat = 0.0
IF (ALLOCATED(jac_rateFactor)) THEN
  DEALLOCATE(jac_rateFactor)
END IF
ALLOCATE(jac_rateFactor(ncomp+nexchange+nsurf))
jac_rateFactor = 0.0
IF (ALLOCATED(jac_pre)) THEN
  DEALLOCATE(jac_pre)
END IF
ALLOCATE(jac_pre(ncomp+nexchange+nsurf,nreactmax))
jac_pre = 0.0
IF (ALLOCATED(jac_preKin)) THEN
  DEALLOCATE(jac_preKin)
END IF
ALLOCATE(jac_preKin(ncomp,nreactkinmax))
jac_preKin = 0.0

IF (ALLOCATED(sppTMP)) THEN
  DEALLOCATE(sppTMP)
END IF
ALLOCATE(sppTMP(ncomp+nspec))
sppTMP = 0.0
IF (ALLOCATED(sppTMP10)) THEN
  DEALLOCATE(sppTMP10)
END IF
ALLOCATE(sppTMP10(ncomp+nspec))
sppTMP10 = 0.0

IF (ALLOCATED(sppTMPperturb)) THEN
  DEALLOCATE(sppTMPperturb)
END IF
ALLOCATE(sppTMPperturb(ncomp+nspec))
sppTMPperturb = 0.0
IF (ALLOCATED(sppTMP10perturb)) THEN
  DEALLOCATE(sppTMP10perturb)
END IF
ALLOCATE(sppTMP10perturb(ncomp+nspec))
sppTMP10perturb = 0.0

!!  Allocation of Crystal Size Distribution stuff

!!CALL AllocateCSD(ncomp,nrct,nx,ny,nz)

!!IF (JennyRifle) THEN
!!
!!  IF (ALLOCATED(tauZero)) THEN
!!    DEALLOCATE(tauZero)
!!  END IF
!!  ALLOCATE(tauZero(nx,ny,nz))
!!  tauZero = 0.0d0
!!
!!  IF (ALLOCATED(MetabolicLag)) THEN
!!    DEALLOCATE(MetabolicLag)
!!  END IF
!!  ALLOCATE(MetabolicLag(nx,ny,nz))
!!  MetabolicLag = 0.0d0
!!
!!END IF

RETURN
END SUBROUTINE AllocateALL