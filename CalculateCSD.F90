SUBROUTINE CalculateCSD(jx,jy,jz,nrct,ncomp,delt)
USE crunchtype
USE params
USE temperature
USE mineral
USE crunch_interface
USE NanoCrystal

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: ncomp
REAL(DP), INTENT(IN)                          :: delt               

! Internal arrays and variables

REAL(DP)                                      :: dxCSD
REAL(DP)                                      :: rinv
REAL(DP)                                      :: fe
REAL(DP)                                      :: fw
REAL(DP)                                      :: ae
REAL(DP)                                      :: aw
REAL(DP)                                      :: apx
REAL(DP)                                      :: pi
REAL(DP)                                      :: Tk
REAL(DP)                                      :: VolumeSingleMolecule
REAL(DP)                                      :: gCritical
REAL(DP)                                      :: qNucleation
REAL(DP)                                      :: qNucleationPrefactor
REAL(DP)                                      :: CriticalRadius
REAL(DP)                                      :: CriticalRadiusCheck
REAL(DP)                                      :: RatioVolumeBoltzmann
REAL(DP)                                      :: dt
REAL(DP),PARAMETER                                      :: AvogadroNumber=6.02214E23
REAL(DP),PARAMETER                                      :: BoltzmannConstant=1.36E-23

INTEGER(I4B)                                  :: l
INTEGER(I4B)                                  :: k
INTEGER(I4B)                                  :: nn


!! Crystal Size Distribution calculation
    
    dt = delt/1.0
    rinv = 1.0/dt
    Tk = t(jx,jy,jz) + 273.15

!!    LinearGrowthRate = 0.0d0

    pi = DACOS(-1.0d0)
    qNucleationPrefactor = 1.0E05
    sigma = 0.0d0
    sigma(3) = 200.0/1000.0
    sigma(4) = 90.0/1000.0

!!  Or calculate the volume of a single molecule from the molar volume divided by Avogradro's number
!! (then multiply by number of units in a cluster)
    
    dxCSD = 2.0E-06
    radius(1) = 4.0E-06
    DO l = 2,nCSD
      radius(l) = radius(l-1) + dxCSD
    END DO

    DO k=1,nrct
 
      IF (CrystalSizeDistribution(k)) THEN

!!       nCrystal(20,k,jx,jy,jz) = 100

!!        LinearGrowthRate = 1.D-8

        VolumeSingleMolecule = volmol(k)/AvogadroNumber
        CALL satcalc(ncomp,nrct,jx,jy,jz)
        RatioVolumeBoltzmann = VolumeSingleMolecule/BoltzmannConstant
        gCritical = (pi*sigma(k)**3 * RatioVolumeBoltzmann*RatioVolumeBoltzmann )/          &
                ( 3.0d0* (Tk*silog(1,k)*clg)**2 )
        qNucleation = qNucleationPrefactor*DEXP (-gCritical/(BoltzmannConstant*Tk) )
        NucleationScaleFactor = qNucleationPrefactor*100000.0
        NucleationScaleFactor = 1.0d0
        qNucleation = qNucleation/NucleationScaleFactor

        DO l = 1,nCSD
          xtvd(l) = DFLOAT(nCrystal(l,k,jx,jy,jz))
        END DO

        fw = LinearGrowthRate(1,k,jx,jy,jz)
        fe = 0.5d0*(LinearGrowthRate(2,k,jx,jy,jz) + LinearGrowthRate(1,k,jx,jy,jz))
        aw = DMAX1(fw,0.0D0)
        ae = DMAX1(-fe,0.0D0)
        apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
        aaCSD(1) = -aw*dt
        bbCSD(1) = apx*dt + dxCSD
        ccCSD(1) = -ae*dt
        uuCSD(1) = xtvd(1)
        rrCSD(1) = xtvd(1)*dxCSD + dt*qNucleation*dxCSD

        fw = 0.5d0*(LinearGrowthRate(nCSD-1,k,jx,jy,jz)+LinearGrowthRate(nCSD,k,jx,jy,jz))
        fe = LinearGrowthRate(nCSD,k,jx,jy,jz)
        aw = DMAX1(fw,0.0D0)
        ae = DMAX1(-fe,0.0D0)
        apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
        aaCSD(nCSD) = -aw*dt
        bbCSD(nCSD) = apx*dt + dxCSD
        ccCSD(nCSD) = -ae*dt
        uuCSD(nCSD) = xtvd(nCSD)
        rrCSD(nCSD) = xtvd(nCSD)*dxCSD

        DO l = 2,nCSD-1

!!        Define fluxes at boundaries
          fw = 0.5d0*(LinearGrowthRate(l-1,k,jx,jy,jz)+LinearGrowthRate(l,k,jx,jy,jz))
          fe = 0.5d0*(LinearGrowthRate(l+1,k,jx,jy,jz)+LinearGrowthRate(l,k,jx,jy,jz))
          aw = DMAX1(fw,0.0D0)
          ae = DMAX1(-fe,0.0D0)
          apx = DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)

          aaCSD(l) = -aw*dt
          bbCSD(l) = apx*dt + dxCSD
          ccCSD(l) = -ae*dt
          uuCSD(l) = xtvd(l)
          rrCSD(l) = xtvd(l)*dxCSD

        END DO
        
        CALL tridag_ser(aaCSD,bbCSD,ccCSD,rrCSD,uuCSD)      
        
        DO l = 1,nCSD
          nCrystal(l,k,jx,jy,jz) = NINT(uuCSD(l))
        END DO

      END IF
          
    END DO   ! End of loop for CSD calculation


RETURN
END SUBROUTINE CalculateCSD