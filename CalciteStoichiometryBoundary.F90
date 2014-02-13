SUBROUTINE CalciteStoichiometryBoundary(nbnd)
USE crunchtype
USE mineral
USE concentration

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                   :: nbnd


REAL(DP)                                                   :: TotUranium
REAL(DP)                                                   :: IsotopeRatio
REAL(DP)                                                   :: Ratio
REAL(DP)                                                   :: tot
REAL(DP)                                                   :: MuTotalUranium
REAL(DP)                                                   :: TotU_Calcite
REAL(DP)                                                   :: TotCa_Calcite

!!  Called at the beginning of the time step to set the calcite stoichiometry for that time step

IsotopeRatio = sbnd(ik234U,nbnd)/(sbnd(ik238U,nbnd)+sbnd(ik234U,nbnd))  !  Total concen.
TotUranium = sbnd(ik234U,nbnd) + sbnd(ik238U,nbnd)
Ratio = DistributionCalcite*TotUranium/sbnd(ikCa,nbnd)                         !  Ratio is mole fraction Uranium relative to Calcium

!!  Ratio is the total moles U in calcite divided by total moles Ca in calcite

TotU_Calcite = 0.97/(1.0D0 + 1.0d0/Ratio)
TotCa_Calcite = TotU_Calcite/Ratio
muTotalUranium = TotU_Calcite/(TotU_Calcite + TotCa_Calcite + 0.03)
muCalciumBoundary(nbnd) = 0.97d0 - muTotalUranium
muUranium234Boundary(nbnd) = IsotopeRatio*muTotalUranium
muUranium238Boundary(nbnd) = (1.0d0 - IsotopeRatio)*muTotalUranium

END SUBROUTINE CalciteStoichiometryBoundary
