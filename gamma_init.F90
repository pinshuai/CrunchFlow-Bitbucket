SUBROUTINE gamma_init(ncomp,nspec,tempc,sqrt_sion)
USE crunchtype
USE params
USE runtime, ONLY: Benchmark
USE concentration
USE temperature

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nspec

REAL(DP), INTENT(IN)                                       :: tempc
REAL(DP), INTENT(OUT)                                      :: sqrt_sion

!  Internal variables

REAL(DP)                                                   :: TotalMoles
REAL(DP)                                                   :: ah
REAL(DP)                                                   :: bh
REAL(DP)                                                   :: bdt
REAL(DP)                                                   :: Chargesum
REAL(DP)                                                   :: sion_tmp
REAL(DP)                                                   :: aa1
REAL(DP)                                                   :: GamWaterCheck

REAL(DP)                                                   :: dhad,dhbd,tempk,tconv

INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: it
INTEGER(I4B)                                               :: ItPoint

ChargeSum = 0.0d0
TotalMoles = 0.0d0
DO ik = 1,ncomp+nspec
  IF (ulab(ik) /= 'H2O') THEN
    TotalMoles = TotalMoles + sptmp10(ik)
    ChargeSum = ChargeSum + sptmp10(ik)*chg(ik)*chg(ik)
  ELSE
    CONTINUE
  END IF
END DO

sion_tmp = 0.50D0*ChargeSum

IF (sion_tmp < 25.0d0) THEN
  sqrt_sion = DSQRT(sion_tmp)
ELSE
  sion_tmp = 0.0d0
  sqrt_sion = 0.0d0
END IF

IF (ntemp == 1) THEN
  ah = adh(1)
  bh = bdh(1)
  bdt = bdot(1)
ELSE

  ItPoint = 0
  DO it = 1,ntemp
    IF (tempc == DatabaseTemperature(it)) THEN
      ItPoint = it
      ah = adh(ItPoint)
      bh = bdh(ItPoint)
      bdt = bdot(ItPoint)
    END IF
  END DO

  IF (ItPoint == 0) THEN
    ah = adhcoeff(1) + adhcoeff(2)*tempc  &
       + adhcoeff(3)*tempc*tempc + adhcoeff(4)*tempc*tempc*tempc  &
       + adhcoeff(5)*tempc*tempc*tempc*tempc
    bh = bdhcoeff(1) + bdhcoeff(2)*tempc  &
       + bdhcoeff(3)*tempc*tempc + bdhcoeff(4)*tempc*tempc*tempc  &
       + bdhcoeff(5)*tempc*tempc*tempc*tempc
    bdt = bdtcoeff(1) + bdtcoeff(2)*tempc  &
       + bdtcoeff(3)*tempc*tempc + bdtcoeff(4)*tempc*tempc*tempc  &
       + bdtcoeff(5)*tempc*tempc*tempc*tempc

  ELSE
    CONTINUE
  END IF

END IF

tconv = 273.15d0
tempk = tempc + tconv
call dhconst(ah,bh,tempk,tconv)

DO ik = 1,ncomp+nspec

  IF (chg(ik) == 0.0d0) THEN

    IF (ulab(ik) == 'H2O') THEN

      gamWaterCheck = 1.0d0 - 0.017d0*TotalMoles
      gamtmp(ik) = DLOG(gamWaterCheck)
    ELSE
        aa1 = 0.10d0*sion_tmp
        gamtmp(IK) = clg*aa1
    END IF

  ELSE

    IF (IncludeBdot) THEN  

      IF (acmp(ik) == 0.0d0) THEN      !! Davies equation
        aa1 = -0.5115d0*chg(IK)*chg(IK)* (sqrt_sion/(1.0d0 + sqrt_sion)  - 0.24d0*sion_tmp  )
      ELSE                  !! WATEQ extended Debye-Huckel

        aa1 = -(ah*chg(IK)*chg(IK)*sqrt_sion)/            &
              (1.0d0 + acmp(IK)*bh*sqrt_sion)                   &         
              + bdotParameter(ik)*sion_tmp
      END IF
    ELSE                    !!  Helgesonian-LLNL bdot expression based on extended Debye-Huckel

      aa1 = -(ah*chg(IK)*chg(IK)*sqrt_sion)/              &
              (1.0d0 + acmp(IK)*bh*sqrt_sion)                  &         
              + bdt*sion_tmp
    END IF

    gamtmp(IK) = clg*aa1

  END IF

END DO

RETURN
END SUBROUTINE gamma_init
!*****************************************************************
