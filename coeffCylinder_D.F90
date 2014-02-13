!  *************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:53:17

SUBROUTINE coeffCylinder_D(nx,ny,nz,ncomp,nspec)
USE crunchtype
USE medium
USE concentration
USE transport
USE temperature
USE params
USE CrunchFunctions

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: nx
INTEGER(I4B), INTENT(IN)                      :: ny
INTEGER(I4B), INTENT(IN)                      :: nz
INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec

! Internal arrays and variables

REAL(DP)                                      :: sum
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
REAL(DP)                                      :: dumsum
REAL(DP)                                      :: avgro
REAL(DP)                                      :: fe
REAL(DP)                                      :: ae
REAL(DP)                                      :: fw
REAL(DP)                                      :: aw
REAL(DP)                                      :: fn
REAL(DP)                                      :: an
REAL(DP)                                      :: fs
REAL(DP)                                      :: as
REAL(DP)                                      :: tk
REAL(DP)                                      :: quirk
REAL(DP)                                      :: porn
REAL(DP)                                      :: pors
REAL(DP)                                      :: dharm
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
REAL(DP)                                      :: zero
REAL(DP)                                      :: de_d
REAL(DP)                                      :: dw_d
REAL(DP)                                      :: ae_d
REAL(DP)                                      :: aw_d
REAL(DP)                                      :: apx_d
REAL(DP)                                      :: dn_d
REAL(DP)                                      :: ds_d
REAL(DP)                                      :: an_d
REAL(DP)                                      :: as_d
REAL(DP)                                      :: apy_d
REAL(DP)                                      :: tort
REAL(DP)                                      :: pi
REAL(DP)                                      :: AreaE
REAL(DP)                                      :: AreaW
REAL(DP)                                      :: AreaS
REAL(DP)                                      :: AreaN

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j
INTEGER(I4B)                                  :: ik

jz = 1

pi = DACOS(-1.0d0)
quirk = 1.0
zero = 0.0D0

IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

DO ik = 1,ncomp+nspec
  d_correct(ik) = d_sp(ik)/d_25
END DO

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    j = (jy-1)*nx+jx
    tk = 273.15 + t(jx,jy,jz)
    IF (idiffus == 0) THEN
      dstar(jx,jy,jz) = dzero*EXP((activation/rgas)*(tk25 - 1.0/tk))/formation
    ELSE
      dstar(jx,jy,jz) = dcoeff/formation
    END IF
  END DO
  dstar(0,jy,jz) = dstar(1,jy,jz)
  dstar(nx+1,jy,jz) = dstar(nx,jy,jz)
END DO
DO jx = 1,nx
  dstar(jx,0,jz) = dstar(jx,1,jz)
  dstar(jx,ny+1,jz) = dstar(jx,ny,jz)
END DO

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    
    porp = por(jx,jy,jz)
    satp = satliq(jx,jy,jz)
    
    IF (nx == 1) GO TO 100  
    
    IF (jx == 1) THEN
      
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*dxx(1)
      pore = por(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (pore > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dume = ro(jx+1,jy,jz)*dstar(jx+1,jy,jz)*sate*pore*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*sate*porp*tort
        dumw = dumpx
      ELSE IF (MillingtonQuirk) THEN
        dume = ro(jx+1,jy,jz)*(sate)**(quirk)*(pore)**(uli)*dstar(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(quirk)*(porp)**(uli)*dstar(jx,jy,jz)
        dumw = dumpx
      ELSE
        dume = ro(jx+1,jy,jz)*sate*pore*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumw = dumpx
      END IF
    ELSE IF (jx == nx) THEN
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      dxe = 0.5d0*dxx(nx)
      pore = por(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx,jy,jz)
      satw = satliq(jx-1,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (porw > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumw = ro(jx-1,jy,jz)*dstar(jx-1,jy,jz)*satw*porw*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort  
        dume = dumpx 
      ELSE IF (MillingtonQuirk) THEN
        dumw = ro(jx-1,jy,jz)*(satw)**(quirk)*(porw)**(uli)*dstar(jx-1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(quirk)*(porp)**(uli)*dstar(jx,jy,jz)
        dume = dumpx
      ELSE
        dumw = ro(jx-1,jy,jz)*satw*porw*dstar(jx-1,jy,jz)*tortuosity(jx-1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dume = dumpx
      END IF
    ELSE
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx-1,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (pore > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dume = ro(jx+1,jy,jz)*dstar(jx+1,jy,jz)*sate*pore*tort
        IF (porw > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumw = ro(jx-1,jy,jz)*dstar(jx-1,jy,jz)*sate*porw*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort
      ELSE IF (MillingtonQuirk) THEN
        dume = ro(jx+1,jy,jz)*(sate)**(quirk)*(pore)**(uli)*dstar(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(quirk)*(porp)**(uli)*dstar(jx,jy,jz)
        dumw = ro(jx-1,jy,jz)*(satw)**(quirk)*(porw)**(uli)*dstar(jx-1,jy,jz)
      ELSE
        dume = ro(jx+1,jy,jz)*sate*pore*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumw = ro(jx-1,jy,jz)*satw*porw*dstar(jx-1,jy,jz)*tortuosity(jx-1,jy,jz)
      END IF
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy == 1) THEN
      dyn = 0.5d0*(dyy(jy)+dyy(jy+1))
      dys = 0.5d0*dyy(1)
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (porn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumn = ro(jx,jy+1,jz)*dstar(jx,jy+1,jz)*satn*porn*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
        dums = dumpy
      ELSE IF (MillingtonQuirk) THEN
        dumn = ro(jx,jy+1,jz)*(satn)**(quirk)*(por(jx,jy+1,jz))**(uli)*dstar(jx,jy+1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(quirk)*(por(jx,jy,jz))**(uli)*dstar(jx,jy,jz)*anisotropyY
        dums = dumpy
      ELSE
        dumn = ro(jx,jy+1,jz)*satn*porn*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
        dumpy = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
        dums = dumpy
      END IF
    ELSE IF (jy == ny) THEN
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      dyn = 0.5d0*dyy(ny)
      porn = por(jx,jy,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy,jz)
      sats = satliq(jx,jy-1,jz)
      IF (UseThresholdPorosity) THEN
        IF (pors > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dums = ro(jx,jy-1,jz)*dstar(jx,jy-1,jz)*sats*pors*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
        dumn = dumpy
      ELSE IF (MillingtonQuirk) THEN
        dums = ro(jx,jy-1,jz)*(sats)**(quirk)*(pors)**(uli)*dstar(jx,jy-1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(quirk)*(por(jx,jy,jz))**(uli)*dstar(jx,jy,jz)*anisotropyY
        dumn = dumpy
      ELSE
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz)*satp*porP*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
        dumn = dumpy
      END IF
    ELSE
      dyn = 0.5d0*(dyy(jy)+dyy(jy+1))
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy-1,jz)

      IF (UseThresholdPorosity) THEN
        IF (porn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumn = ro(jx,jy+1,jz)*dstar(jx,jy+1,jz)*satn*porn*tort*anisotropyY
        IF (pors > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dums = ro(jx,jy-1,jz)*dstar(jx,jy-1,jz)*sats*pors*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
      ELSE IF (MillingtonQuirk) THEN
        dumn = ro(jx,jy+1,jz)*(satn)**(quirk)*(porn)**(uli)*dstar(jx,jy+1,jz)*anisotropyY
        dums = ro(jx,jy-1,jz)*(sats)**(quirk)*(pors)**(uli)*dstar(jx,jy-1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(quirk)*(porp)**(uli)*dstar(jx,jy,jz)*anisotropyY
      ELSE
        dumn = ro(jx,jy+1,jz)*satn*porn*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
      END IF
    END IF
    
    200     CONTINUE
    IF (nx == 1) GO TO 300
    
    IF (jx == 1) THEN
      
      avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
      
 !!     dharm = ArithmeticMean(dume,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      dspe = avgro*dspx(jx,jy,jz)
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      ae = DMAX1(-fe,zero) + de
      netflowx(1,jy,jz) = qx(jx,jy,jz) + FluidBuryX(jx)
      de_d = AreaE*dharm/dxe
      ae_d = de_d
      
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)
      avgro = ro(jx,jy,jz)
      dharm = dumpx
      dspw = avgro*dspx(jx,jy,jz)
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(jx-1,jy,jz) + FluidBuryX(jx-1))
      IF (jc(1) == 2) THEN
        aw = DMAX1(fw,zero)      ! Pure advective boundary
      ELSE
        aw = DMAX1(fw,zero) + dw
      END IF
      IF (jc(1) == 2) THEN
        dw_d = 0.0
        aw_d = 0.0
      ELSE
        dw_d = AreaW*dharm/dxw
        aw_d = dw_d
      END IF
      
      netflowx(0,jy,jz) = qx(jx-1,jy,jz) + FluidBuryX(jx-1)
      
      IF (jc(1) == 2) THEN
        apx = de + DMAX1(-fw,zero) + DMAX1(fe,zero)      !  Pure advective boundary
        apx_d = de_d 
      ELSE
        apx = dw + de + DMAX1(-fw,zero) + DMAX1(fe,zero)
        apx_d = de_d + dw_d
      END IF

      
    ELSE IF (jx == nx) THEN
      
      avgro = ro(jx,jy,jz)
      dharm = dumpx
      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      dspe = avgro*dspx(jx,jy,jz)
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      IF (jc(2) == 2) THEN
        ae = DMAX1(-fe,zero)
      ELSE
        ae = DMAX1(-fe,zero) + de
      END IF
      netflowx(jx,jy,jz) = qx(jx,jy,jz) + FluidBuryX(jx)
      IF (jc(2) == 2) THEN
        de_d = 0.0
        ae_d = de_d
      ELSE
        de_d = AreaE*dharm/dxe
        ae_d = de_d
      END IF
      
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dumw,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)
      dspw = avgro*dspx(jx-1,jy,jz)
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(jx-1,jy,jz) + FluidBuryX(jx-1))
      aw = DMAX1(fw,zero) + dw
      dw_d = AreaW*dharm/dxw
      aw_d = dw_d
      
      IF (jc(2) == 2) THEN
        apx = dw + DMAX1(-fw,zero) + DMAX1(fe,zero)
        apx_d = dw_d
      ELSE
        apx = dw + de + DMAX1(-fw,zero) + DMAX1(fe,zero)
        apx_d = de_d + dw_d
      END IF
      
    ELSE
      
      avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dume,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      dspe = avgro*dspx(jx,jy,jz)
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      ae = DMAX1(-fe,zero) + de
      netflowx(1,jy,jz) = qx(jx,jy,jz) + FluidBuryX(jx)
      de_d = AreaE*dharm/dxe
      ae_d = de_d
      
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dumw,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)
      dspw = avgro*dspx(jx-1,jy,jz)
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(jx-1,jy,jz) + FluidBuryX(jx-1))
      aw = DMAX1(fw,zero) + dw
      dw_d = AreaW*dharm/dxw
      aw_d = dw_d
      
      apx = dw + de + DMAX1(-fw,zero) + DMAX1(fe,zero)
      apx_d = de_d + dw_d
      
    END IF
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    IF (jy == 1) THEN
      
      avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dumn,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )         !! Area to integrate qx over (this assumes that X is the radial direction):  pi*r^2_(jx+1) - pi*r^2_(jx-1)
      dspn = avgro*dspy(jx,jy,jz)
      dn = AreaN*dspn/dyn                                                        !! Integrated dispersive flux
      fn = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))                           !! Integrated advective flux
      an = DMAX1(-fn,zero) + dn
      dn_d = AreaN*dharm/dyn                                                     !! Integrated diffusive flux
      an_d = dn_d
      
      avgro = ro(jx,jy,jz)
      dharm = dumpy
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      dsps = avgro*dspy(jx,jy,jz) 
      ds = AreaS*dsps/dys
      fs = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))

      IF (jc(3) == 2) THEN
        as = DMAX1(fs,zero)
      ELSE
        as = DMAX1(fs,zero) + ds
      END IF
      
      IF (jc(3) == 2) THEN
        ds_d = 0.0
        as_d = 0.0
      ELSE
        ds_d = AreaS*dharm/dys
        as_d = ds_d
      END IF

      IF (jc(3) == 2) THEN
        apy = dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
        apy_d = dn_d 
      ELSE
        apy = ds + dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
        apy_d = dn_d + ds_d
      END IF
      
    ELSE IF (jy == ny) THEN
      
      avgro = ro(jx,jy,jz)
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      dspn = avgro*dspy(jx,jy-1,jz) 
      dn = AreaN*dspn/dyn
      fn = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))
      IF (jc(4) == 2) THEN
        an = DMAX1(-fn,zero)
      ELSE
        an = DMAX1(-fn,zero) + dn
      END IF
      dharm = dumpy
      IF (jc(4) == 2) THEN
        dn_d = 0.0
        an_d = 0.0
      ELSE
        dn_d = AreaN*dharm/dyn
        an_d = dn_d
      END IF
      
      avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dums,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      dsps = avgro*dspy(jx,jy-1,jz) 
      ds = AreaS*dsps/dys
      fs = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))
      as = DMAX1(fs,zero) + ds
      ds_d = AreaS*dharm/dys
      as_d = ds_d

      IF (jc(4) == 2) THEN
        apy = ds + DMAX1(-fs,zero) + DMAX1(fn,zero)
        apy_d = ds_d 
      ELSE
        apy = ds + dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
        apy_d = dn_d + ds_d
      END IF
      
    ELSE
      
      avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
!!      dharm = ArithmeticMean(dumn,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      dspn = avgro*dspy(jx,jy,jz)
      dn = AreaN*dspn/dyn
      fn = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))
      an = DMAX1(-fn,zero) + dn
      dn_d = AreaN*dharm/dyn
      an_d = dn_d
      
      avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
!!       dharm = ArithmeticMean(dums,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      dsps = avgro*dspy(jx,jy-1,jz) 
      ds = AreaS*dsps/dys
      fs = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))
      as = DMAX1(fs,zero) + ds
      ds_d = AreaS*dharm/dys
      as_d = ds_d
      
      apy = ds + dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
      apy_d = dn_d + ds_d
      
    END IF
    
    400     CONTINUE
    IF (nx == 1) THEN
      a(jx,jy,jz) = zero
      b(jx,jy,jz) = zero
      c(jx,jy,jz) = zero
    ELSE
      
      a(jx,jy,jz) = -aw
      c(jx,jy,jz) = -ae
      b(jx,jy,jz) = apx
      a_d(jx,jy,jz) = -aw_d
      c_d(jx,jy,jz) = -ae_d
      b_d(jx,jy,jz) = apx_d
      
    END IF
    
    IF (ny == 1) THEN
      d(jx,jy,jz) = zero
      f(jx,jy,jz) = zero
      e(jx,jy,jz) = zero
    ELSE
      
      d(jx,jy,jz) = -an
      f(jx,jy,jz) = -as
      e(jx,jy,jz) = apy
      d_d(jx,jy,jz) = -an_d
      f_d(jx,jy,jz) = -as_d
      e_d(jx,jy,jz) = apy_d
      
    END IF
    dxy(jx,jy,jz) = pi*dyy(jy)*( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )          !! Volume in 2D system with c cylindrical coordinates (assumes X is radial direction)
    
  END DO
END DO


RETURN
END SUBROUTINE coeffCylinder_D
