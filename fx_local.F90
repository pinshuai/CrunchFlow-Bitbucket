!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:43
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE fx_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,  &
    nspec,ngas,neqn,dt,jx,jy,jz,nx,ny,nz,ne,time)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature
USE RunTime

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

integer nx,ny,nz,ne
INTEGER(I4B), INTENT(IN)                                  :: ncomp
INTEGER(I4B), INTENT(IN)                                  :: nexchange
INTEGER(I4B), INTENT(IN)                                  :: nexch_sec
INTEGER(I4B), INTENT(IN)                                  :: nsurf
INTEGER(I4B), INTENT(IN)                                  :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                  :: nrct
INTEGER(I4B), INTENT(IN)                                  :: nspec
INTEGER(I4B), INTENT(IN)                                  :: ngas
INTEGER(I4B), INTENT(IN)                                  :: neqn
INTEGER(I4B), INTENT(IN)                                  :: jx
INTEGER(I4B), INTENT(IN)                                  :: jy
INTEGER(I4B), INTENT(IN)                                  :: jz

REAL(DP), INTENT(IN)                                      :: dt
REAL(DP), INTENT(IN)                                      :: time

!  Internal variables and arrays

REAL(DP), DIMENSION(ncomp)                                :: scorr
REAL(DP)                                                  :: r
REAL(DP)                                                  :: satgasnew
REAL(DP)                                                  :: satgasold
REAL(DP)                                                  :: df
REAL(DP)                                                  :: source
REAL(DP)                                                  :: aq_accum
REAL(DP)                                                  :: gas_accum
REAL(DP)                                                  :: ex_accum
REAL(DP)                                                  :: surf_accum
REAL(DP)                                                  :: satl
REAL(DP)                                                  :: satlold
REAL(DP)                                                  :: Retardation

INTEGER(I4B)                                              :: i
INTEGER(I4B)                                              :: is
INTEGER(I4B)                                              :: ix
INTEGER(I4B)                                              :: ind

REAL(DP)                                                  :: CellVolume
REAL(DP)                                                  :: MultiplyCell
REAL(DP)                                                  :: pi

!  This routine calculates (and assembles) the function residuals.

pi = DACOS(-1.0d0)

r = 1.0D0/dt
IF (cylindrical) THEN
  df = 1.0/(3.1416*dxx(jx)*dxx(jx)*dzz(jx,jy,jz))
  CellVolume = dyy(jy)*pi*( (x(jx)+dxx(jx)/2.0d0 )**2.0d0 - ( x(jx)-dxx(jx)/2.0d0 )**2.0d0  )
  df = 1.0d0
  MultiplyCell = CellVolume
  IF (CylindricalDivideVolume) THEN
      df = 1.0/CellVolume
      MultiplyCell = 1.0
  END IF
ELSE
  df = 1.0/(dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
END IF

!!   CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
!!   MultiplyCell = CellVolume

Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

fxmax = 0.0
satl = satliq(jx,jy,jz)
satgasnew = 1.0 - satl
satlOld = satl
!!satlold = satliqold(jx,jy,jz)
satgasold = 1.0 - satlold

DO i = 1,ncomp  

  ind = i

  source = 0.0
  
!!      Assumes Kd is dimensionless  

!!  aq_accum = xgram(jx,jy,jz)*r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
!!           (satl*s(i,jx,jy,jz)-satlold*sn(i,jx,jy,jz)) *(1.0 + Retardation*distrib(i) ) 

  if (i /= ikh2o) THEN
    aq_accum = r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
           (H2Oreacted(jx,jy,jz)*xgram(jx,jy,jz)*satl*s(i,jx,jy,jz) - &
            xgramOld(jx,jy,jz)*satlold*sn(i,jx,jy,jz)) *(1.0 + Retardation*distrib(i) ) 
    
  ELSE
    aq_accum = r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
           (xgram(jx,jy,jz)*satl*s(i,jx,jy,jz) - &
            xgramOld(jx,jy,jz)*satlold*sn(i,jx,jy,jz)) *(1.0 + Retardation*distrib(i) ) 
  END IF

  IF (isaturate == 1) THEN
    gas_accum = por(jx,jy,jz)*r*(satgasnew*sgas(i,jx,jy,jz)-satgasold*sgasn(i,jx,jy,jz))
  ELSE
    gas_accum = 0.0
  END IF

  ex_accum = r*(sNCexch_local(i) - sexold(i,jx,jy,jz))      &
            +    r*(sNCsurf_local(i) - ssurfold(i,jx,jy,jz))

  fxx(ind) = aq_accum + gas_accum + ex_accum - source
  
  IF (ABS(fxx(ind)) > fxmax(i)) THEN
    fxmax(i) = ABS(fxx(ind))
  END IF
  
END DO

DO is = 1,nsurf
  ind = is+ncomp+nexchange
  surf_accum = r*(ssurf_local(is) - ssurfn(is,jx,jy,jz))
  fxx(ind) = surf_accum
END DO

RETURN
END SUBROUTINE fx_local
!***************************************************************************
