SUBROUTINE CalciteStoichiometry(nx,ny,nz)
USE crunchtype
USE mineral
USE concentration

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

REAL(DP)                                                   :: TotUranium
REAL(DP)                                                   :: IsotopeRatio
REAL(DP)                                                   :: Ratio
REAL(DP)                                                   :: tot
REAL(DP)                                                   :: MuTotalUranium
REAL(DP)                                                   :: TotU_Calcite
REAL(DP)                                                   :: TotCa_Calcite

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

!!  Called at the beginning of the time step to set the calcite stoichiometry for that time step

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      IsotopeRatio = s(ik234U,jx,jy,jz)/(s(ik238U,jx,jy,jz)+s(ik234U,jx,jy,jz))  !  Total concen.

      TotUranium = s(ik234U,jx,jy,jz) + s(ik238U,jx,jy,jz)

      Ratio = DistributionCalcite*TotUranium/s(ikCa,jx,jy,jz)                           !  Ratio is mole fraction Uranium relative to Calcium

!!  Ratio is the total moles U in calcite divided by total moles Ca in calcite

      TotU_Calcite = 0.97/(1.0D0 + 1.0d0/Ratio)
      TotCa_Calcite = TotU_Calcite/Ratio
      muTotalUranium = TotU_Calcite/(TotU_Calcite + TotCa_Calcite + 0.03)
      muCalcium(jx,jy,jz) = 0.97d0 - muTotalUranium
      muUranium234(jx,jy,jz) = IsotopeRatio*muTotalUranium
      muUranium238(jx,jy,jz) = (1.0d0 - IsotopeRatio)*muTotalUranium

    END DO
  END DO
END DO


END SUBROUTINE CalciteStoichiometry
