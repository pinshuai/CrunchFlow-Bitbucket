SUBROUTINE ExchangeConcentration(jx,jy,jz,ncomp,nexchange,nexch_sec)
USE crunchtype
USE concentration
USE RunTime
USE temperature
USE transport
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nexchange
INTEGER(I4B), INTENT(IN)                                    :: nexch_sec

!  Internal variables and arrays

REAL(DP)                                                    :: AqueousToBulk 
INTEGER(I4B)                                                :: i,ns

!! NOTE:  Assumes a "bare" primary species site (cycling only "secondary" surface complexes)

CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

ExchangeCon = 0.0

DO i = 1,ncomp  
  DO ns = 1,nexch_sec
    ExchangeCon(i) = ExchangeCon(i) + muexc(ns,i)*spex10(ns+nexchange,jx,jy,jz)
  END DO

  IF (AqueousToBulk /= 0.0) THEN        !!  Completely dry cell
    ExchangeCon(i) = ExchangeCon(i)/AqueousToBulk
  END IF
END DO

RETURN
END SUBROUTINE ExchangeConcentration
