!  *************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:58:01

SUBROUTINE gasdiff(nx,ny,nz)
USE crunchtype
USE params
USE concentration
USE medium
USE transport
USE flow
USE temperature

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: nx
INTEGER(I4B), INTENT(IN)                      :: ny
INTEGER(I4B), INTENT(IN)                      :: nz

! Internal arrays and variables

REAL(DP)                                      :: porp
REAL(DP)                                      :: pore
REAL(DP)                                      :: porw
REAL(DP)                                      :: satw
REAL(DP)                                      :: satp
REAL(DP)                                      :: sate
REAL(DP)                                      :: satn
REAL(DP)                                      :: sats
REAL(DP)                                      :: dumw
REAL(DP)                                      :: dume
REAL(DP)                                      :: dumpx
REAL(DP)                                      :: dumn
REAL(DP)                                      :: dums
REAL(DP)                                      :: dumpy
REAL(DP)                                      :: dxe
REAL(DP)                                      :: dxw
REAL(DP)                                      :: dys
REAL(DP)                                      :: dyn
REAL(DP)                                      :: fe
REAL(DP)                                      :: ae
REAL(DP)                                      :: fw
REAL(DP)                                      :: aw
REAL(DP)                                      :: fn
REAL(DP)                                      :: an
REAL(DP)                                      :: fs
REAL(DP)                                      :: as
REAL(DP)                                      :: tk
REAL(DP)                                      :: QuirkGas
REAL(DP)                                      :: porn
REAL(DP)                                      :: pors
REAL(DP)                                      :: dspe
REAL(DP)                                      :: dspw
REAL(DP)                                      :: de
REAL(DP)                                      :: dw
REAL(DP)                                      :: ds
REAL(DP)                                      :: dn
REAL(DP)                                      :: dsps
REAL(DP)                                      :: dspn
REAL(DP)                                      :: apx
REAL(DP)                                      :: apy
REAL(DP)                                      :: UliGas
REAL(DP)                                      :: zero
REAL(DP)                                      :: gasd
REAL(DP)                                      :: tempe
REAL(DP)                                      :: tempw
REAL(DP)                                      :: temps
REAL(DP)                                      :: tempn

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j

!!IF (MillingtonQuirk) THEN
!!  UliGas = 7.0d0/3.0d0
!!  QuirkGas = 1.0d0/3.0d0
!!ELSE
!!  UliGas = 0.0d0
!!  QuirkGas = 0.0d0
!!END IF

UliGas = 1.0d0
QuirkGas = 0.0d0
zero = 0.0d0

  UliGas = 7.0d0/3.0d0
  QuirkGas = 1.0d0/3.0d0

jz = 1
DO jy = 1,ny
  DO jx = 1,nx

    j = (jy-1)*nx+jx
    tk = 273.15 + t(jx,jy,jz)
    
    porp = por(jx,jy,jz)
    satp = 1.0-satliq(jx,jy,jz)
    
    IF (nx == 1) GO TO 100
    
    IF (jx == 1) THEN
      dxe = 0.5*(dxx(jx)+dxx(jx+1))
      dxw = 0.5*dxx(1)
      pore = por(jx+1,jy,jz)
      sate = 1.0-satliq(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      satw = 1.0-satliq(jx,jy,jz)
      gasd = (pore)**QuirkGas*(sate)**(UliGas)*dgas
      dume = pore*sate*gasd
      gasd = (porp)**QuirkGas*(satp)**(UliGas)*dgas
      dumpx = porp*satp*gasd
      dumw = dumpx
    ELSE IF (jx == nx) THEN
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      dxe = 0.5*dxx(nx)
      pore = por(jx,jy,jz)
      sate = 1.0-satliq(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      satw = 1.0-satliq(jx-1,jy,jz)
      gasd = (porw)**QuirkGas*(satw)**(UliGas)*dgas
      dumw = porw*satw*gasd
      gasd = porp**QuirkGas*(satp)**(UliGas)*dgas
      dumpx = porp*satp*gasd
      dume = dumpx
    ELSE
      dxe = 0.5*(dxx(jx)+dxx(jx+1))
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      sate = 1.0-satliq(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      satw = 1.0-satliq(jx-1,jy,jz)
      gasd = (pore)**QuirkGas*(sate)**(UliGas)*dgas
      dume = pore*sate*gasd
      gasd = (porw)**QuirkGas*(satw)**(UliGas)*dgas
      dumw = porw*satw*gasd
      gasd = (porp)**QuirkGas*(satp)**(UliGas)*dgas
      dumpx = porp*satp*gasd
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy == 1) THEN
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*dyy(1)
      porn = por(jx,jy+1,jz)
      satn = 1.0-satliq(jx,jy+1,jz)
      pors = por(jx,jy,jz)
      sats = 1.0-satliq(jx,jy,jz)
      gasd = (porn)**QuirkGas*(satn)**(UliGas)*dgas
      dumn = porn*satn*gasd
      gasd = (porp)**QuirkGas*(satp)**(UliGas)*dgas
      dumpy = porp*satp*gasd
      dums = dumpy
    ELSE IF (jy == ny) THEN
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      dyn = 0.5*dyy(ny)
      porn = por(jx,jy,jz)
      satn = 1.0-satliq(jx,jy,jz)
      pors = por(jx,jy-1,jz)
      sats = 1.0-satliq(jx,jy-1,jz)
      gasd = (pors)**QuirkGas*(sats)**(UliGas)*dgas
      dums = pors*sats*gasd
      gasd = (porp)**QuirkGas*(satp)**(UliGas)*dgas
      dumpy = porp*satp*gasd
      dumn = dumpy
    ELSE
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      satn = 1.0-satliq(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
      sats = 1.0-satliq(jx,jy-1,jz)
      gasd = (pors)**QuirkGas*(sats)**(UliGas)*dgas
      dums = pors*sats*gasd
      gasd = (porn)**QuirkGas*(satn)**(UliGas)*dgas
      dumn = porn*satn*gasd
      gasd = (porp)**QuirkGas*(satp)**(UliGas)*dgas
      dumpy = porp*satp*gasd
    END IF
    
    200     CONTINUE
    IF (nx == 1) GO TO 300
    
    tempe = dume + dumpx
    tempw = dumw + dumpx
    IF (tempe /= zero) THEN
      dspe = 2.0D0*dume*dumpx/(tempe)
    ELSE
      dspe = zero
    END IF
    IF (tempw /= zero) THEN
      dspw = 2.d0*dumw*dumpx/(tempw)
    ELSE
      dspw = zero
    END IF
    de = dspe*dyy(jy)/dxe
    dw = dspw*dyy(jy)/dxw


    IF (jx == 1 .AND. jc(1) == 2) THEN
      dw = 0.00
    END IF

    IF (jx == NX .AND. jc(2) == 2) THEN
      de = 0.00
    END IF
    
!!  The following forces Dirichlet conditions for gases at all times (unless above is uncommented)

    fe = dyy(jy)*qxgas(jx,jy,jz)
    fw = dyy(jy)*qxgas(jx-1,jy,jz)

    ae = DMAX1(-fe,zero) + de
    aw = DMAX1(fw,zero) + dw
    apx = dw + de + DMAX1(-fw,zero) + DMAX1(fe,zero)
!!  **************************************************************

      netflowx(0,jy,jz) = qx(jx-1,jy,jz) + FluidBuryX(jx-1)
      
      IF (jc(1) == 2) THEN
        apx = de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)  !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    tempn= dumn + dumpy
    temps = dums + dumpy
    IF (tempn /= zero) THEN
      dspn = 2.0D0*dumn*dumpy/(tempn)
    ELSE
      dspn = zero
    END IF
    IF (temps /= zero) THEN
      dsps = 2.0D0*dums*dumpy/(temps)
    ELSE
      dsps = zero
    END IF
    dn = dspn*dxx(jx)/dyn
    ds = dsps*dxx(jx)/dys

    IF (jy == 1 .AND. jc(3) == 2) THEN
      ds = 0.00
    END IF

    IF (jy == Ny .AND. jc(4) == 2) THEN
      dn = 0.00
    END IF
    
!!  The following forces Dirichlet conditions for gases at all times (unless above is uncommented)

    fn = dxx(jx)*qygas(jx,jy,jz)
    an = DMAX1(-fn,zero) + dn
    
    fs = dxx(jx)*qygas(jx,jy-1,jz)
    as = DMAX1(fs,zero) + ds
    
    apy = ds + dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
    
    400     CONTINUE
    
    IF (nx == 1) THEN
      ag(jx,jy,jz) = zero
      cg(jx,jy,jz) = zero
      bg(jx,jy,jz) = zero
    ELSE
      
      ag(jx,jy,jz) = -aw
      cg(jx,jy,jz) = -ae
      bg(jx,jy,jz) = apx
    END IF
    
    IF (ny == 1) THEN
      dg(jx,jy,jz) = zero
      fg(jx,jy,jz) = zero
      eg(jx,jy,jz) = zero
    ELSE
      dg(jx,jy,jz) = -an
      fg(jx,jy,jz) = -as
      eg(jx,jy,jz) = apy
    END IF
    
  END DO
END DO

RETURN
END SUBROUTINE gasdiff
