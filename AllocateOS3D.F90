SUBROUTINE AllocateOS3D(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
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

IMPLICIT NONE

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

  IF (ALLOCATED(cflowin)) THEN
    DEALLOCATE(cflowin)
  END IF
  ALLOCATE(cflowin(ncomp))
  IF (ALLOCATED(cflowout)) THEN
    DEALLOCATE(cflowout)
  END IF
  ALLOCATE(cflowout(ncomp))
  cflowin = 0.0
  cflowout = 0.0
  IF (ALLOCATED(ftvd)) THEN
    DEALLOCATE(ftvd)
  END IF
  ALLOCATE(ftvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(gtvd)) THEN
    DEALLOCATE(gtvd)
  END IF
  ALLOCATE(gtvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(htvd)) THEN
    DEALLOCATE(htvd)
  END IF
  ALLOCATE(htvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(sorp)) THEN
    DEALLOCATE(sorp)
  END IF
  ALLOCATE(sorp(nx,ny,nz))
  sorp = 1.0
  IF (ALLOCATED(gam)) THEN
    DEALLOCATE(gam)
  END IF
  ALLOCATE(gam(ncomp+nspec,nx,ny,nz))
  IF (ALLOCATED(keqaq)) THEN
    DEALLOCATE(keqaq)
  END IF
  ALLOCATE(keqaq(nspec,nx,ny,nz))
  IF (ALLOCATED(keqgas)) THEN
    DEALLOCATE(keqgas)
  END IF
  ALLOCATE(keqgas(ngas,nx,ny,nz))
  IF (ALLOCATED(keqsurf)) THEN
    DEALLOCATE(keqsurf)
  END IF
  ALLOCATE(keqsurf(nsurf_sec,nx,ny,nz))
  IF (ALLOCATED(keqmin)) THEN
    DEALLOCATE(keqmin)
  END IF
  ALLOCATE(keqmin(nreactmax,nrct,nx,ny,nz))
  IF (ALLOCATED(fjac_loc)) THEN
    DEALLOCATE(fjac_loc)
  END IF
  ALLOCATE(fjac_loc(ncomp,ncomp))
  gam = 0.0
  keqaq = 0.0
  keqgas = 0.0
  keqsurf = 0.0
  keqmin = 0.0
  fjac_loc = 0.0
  IF (ALLOCATED(ssurfn)) THEN
    DEALLOCATE(ssurfn)
  END IF
  ALLOCATE(ssurfn(nsurf,nx,ny,nz))
  IF (ALLOCATED(LogTotalSurface)) THEN
    DEALLOCATE(LogTotalSurface)
  END IF
  ALLOCATE(LogTotalSurface(nsurf,nx,ny,nz))
  IF (ALLOCATED(sexold)) THEN
    DEALLOCATE(sexold)
  END IF
  ALLOCATE(sexold(ncomp,nx,ny,nz))
  IF (ALLOCATED(ssurfold)) THEN
    DEALLOCATE(ssurfold)
  END IF
  ALLOCATE(ssurfold(ncomp,nx,ny,nz))
  IF (ALLOCATED(skdold)) THEN
    DEALLOCATE(skdold)
  END IF
  ALLOCATE(skdold(ncomp,nx,ny,nz))
  ssurfn = 0.0
  sexold = 0.0
  ssurfold = 0.0
  skdold = 0.0
  IF (ALLOCATED(fxx)) THEN
    DEALLOCATE(fxx)
  END IF
  ALLOCATE(fxx(neqn))
  fxx = 0.0
  IF (ALLOCATED(sNCexch_local)) THEN
    DEALLOCATE(sNCexch_local)
  END IF
  ALLOCATE(sNCexch_local(ncomp))
  IF (ALLOCATED(sNCsurf_local)) THEN
    DEALLOCATE(sNCsurf_local)
  END IF
  ALLOCATE(sNCsurf_local(ncomp))
  IF (ALLOCATED(ssurf_local)) THEN
    DEALLOCATE(ssurf_local)
  END IF
  ALLOCATE(ssurf_local(nsurf))
  IF (ALLOCATED(fch_local)) THEN
    DEALLOCATE(fch_local)
  END IF
  ALLOCATE(fch_local(ncomp,neqn))
  IF (ALLOCATED(fsurf_local)) THEN
    DEALLOCATE(fsurf_local)
  END IF
  ALLOCATE(fsurf_local(nsurf,nsurf+ncomp))
  ssurf_local = 0.0
  fch_local = 0.0
  fsurf_local = 0.0
  sNCexch_local = 0.0
  sNCsurf_local = 0.0
  IF (ALLOCATED(aaa)) THEN
    DEALLOCATE(aaa)
  END IF
  ALLOCATE(aaa(neqn,neqn))
  IF (ALLOCATED(bb)) THEN
    DEALLOCATE(bb)
  END IF
  ALLOCATE(bb(neqn))
  IF (ALLOCATED(xn)) THEN
    DEALLOCATE(xn)
  END IF
  ALLOCATE(xn(neqn))
  IF (ALLOCATED(indd)) THEN
    DEALLOCATE(indd)
  END IF
  ALLOCATE(indd(neqn))
  aaa = 0.0
  bb = 0.0
  xn = 0.0
  indd = 0
  IF (ALLOCATED(ctvd)) THEN
    DEALLOCATE(ctvd)
  END IF
  ALLOCATE(ctvd(-1:nx+2,-1:ny+2,-1:nz+2))
  ctvd = 0.0
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
  END IF
  ALLOCATE(LogPotential(nsurf,nx,ny,nz))
  LogPotential = 0.0
  IF (ALLOCATED(surfcharge)) THEN
    DEALLOCATE(surfcharge)
  END IF
  ALLOCATE(surfcharge(nrct))
  surfcharge = 0.0
  IF (ALLOCATED(fjpotncomp_local)) THEN
    DEALLOCATE(fjpotncomp_local)
  END IF
  ALLOCATE(fjpotncomp_local(npot,ncomp))
  IF (ALLOCATED(fjpotnsurf_local)) THEN
    DEALLOCATE(fjpotnsurf_local)
  END IF
  ALLOCATE(fjpotnsurf_local(npot,nsurf))
  fjpotncomp_local = 0.0
  fjpotnsurf_local = 0.0
  IF (spherical) THEN
    IF (ALLOCATED(a)) THEN
      DEALLOCATE(a)
    END IF
    ALLOCATE(a(nx,ny,nz))
    IF (ALLOCATED(b)) THEN
      DEALLOCATE(b)
    END IF
    ALLOCATE(b(nx,ny,nz))
    IF (ALLOCATED(c)) THEN
      DEALLOCATE(c)
    END IF
    ALLOCATE(c(nx,ny,nz))
    IF (ALLOCATED(d)) THEN
      DEALLOCATE(d)
    END IF
    ALLOCATE(d(nx,ny,nz))
    IF (ALLOCATED(e)) THEN
      DEALLOCATE(e)
    END IF
    ALLOCATE(e(nx,ny,nz))
    IF (ALLOCATED(f)) THEN
      DEALLOCATE(f)
    END IF
    ALLOCATE(f(nx,ny,nz))
    a = 0.0
    b = 0.0
    c = 0.0
    d = 0.0
    e = 0.0
    f = 0.0
  END IF

!!IF (nsurf > 0) THEN
  IF (ALLOCATED(SurfaceCon)) THEN
    DEALLOCATE(SurfaceCon)
  END IF
  ALLOCATE(SurfaceCon(ncomp))
!!END IF
!!IF (nexchange > 0) THEN
  IF (ALLOCATED(ExchangeCon)) THEN
    DEALLOCATE(ExchangeCon)
  END IF
  ALLOCATE(ExchangeCon(ncomp))
  IF (ALLOCATED(fweight)) THEN
    DEALLOCATE(fweight)
  END IF
  ALLOCATE(fweight(neqn,neqn))
!!END IF

IF (ALLOCATED(IntegerArray)) THEN
  DEALLOCATE(IntegerArray)
END IF
ALLOCATE(IntegerArray(nx,ny,nz))
IF (ALLOCATED(RealArray)) THEN
  DEALLOCATE(RealArray)
END IF
ALLOCATE(RealArray(nx,ny,nz))

RETURN
END SUBROUTINE AllocateOS3D