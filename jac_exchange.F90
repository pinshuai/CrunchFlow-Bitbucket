
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:59:20
 
!******************        GIMRT98     ************************
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

SUBROUTINE jac_exchange(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,neqn,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver
USE transport
USE medium
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                           :: neqn
INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: jx
INTEGER(I4B), INTENT(IN)                           :: jy
INTEGER(I4B), INTENT(IN)                           :: jz

!  Internal variables and arrays

REAL(DP)                                           :: mutemp
REAL(DP)                                           :: spex_conc
REAL(DP)                                           :: sum1
REAL(DP)                                           :: sum2

INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: i2
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ix
INTEGER(I4B)                                       :: ixcheck
INTEGER(I4B)                                       :: ns
INTEGER(I4B)                                       :: is2

fch(:,:,jx,jy,jz) = 0.0

IF (iexc == 1 .OR. iexc == 3) THEN  ! Gapon or Gaines-Thomas (equivalent fractions)
  
  DO nex = 1,nexch_sec
    spex_conc = spex10(nex+nexchange,jx,jy,jz)
    DO i = 1,ncomp
      IF (muexc(nex,i) /= 0.0) THEN
        mutemp = muexc(nex,i) 
        DO i2 = 1,ncomp+nexchange
          fch(i,i2,jx,jy,jz) = fch(i,i2,jx,jy,jz) +    &
                mutemp*muexc(nex,i2)*spex_conc
        END DO
      END IF
    END DO
  END DO

ELSE    !                  Vanselow convention
  
!  fexch = 0.0
  fweight = 0.0
  DO nex = 1,nexch_sec
    spex_conc = aexch(nex)
!    ixcheck = ixlink(nex)
    DO ix = 1,nexchange
      IF (muexc(nex,ix+ncomp) /= 0.0) THEN
        mutemp = muexc(nex,ix+ncomp)
        DO i2 = 1,ncomp+nexchange
          fweight(ix+ncomp,i2) = fweight(ix+ncomp,i2) +     &
              mutemp*muexc(nex,i2)*spex_conc
        END DO
      END IF
!      IF (ixcheck == ix) THEN
!        DO i2 = 1,ncomp+nexchange
!          fexch(ix+ncomp,i2) = fexch(ix+ncomp,i2) +       & 
!              muexc(nex,i2)*spex_conc
!        END DO
!      END IF
    END DO
  END DO
  
  DO nex = 1,nexch_sec
    spex_conc = aexch(nex)
    ix = ixlink(nex)
    DO i = 1,ncomp
      IF (muexc(nex,i) /= 0.0) THEN
        mutemp = muexc(nex,i)
        DO i2 = 1,ncomp+nexchange
          sum1 = mutemp*muexc(nex,i2)*spex_conc*tec(ix)
          sum2 = -exchangesites(ix,jx,jy,jz)*mutemp*spex_conc/  &
             (wt_aexch(ix)*wt_aexch(ix))*fweight(ix+ncomp,i2)
          fch(i,i2,jx,jy,jz) = fch(i,i2,jx,jy,jz) + sum1 + sum2
        END DO
      END IF
    END DO
  END DO
  
END IF

DO ns = 1,nsurf_sec
  spex_conc = spsurf10(ns+nsurf,jx,jy,jz)
  DO i = 1,ncomp
    IF (musurf(ns,i) /= 0.0) THEN
      mutemp = musurf(ns,i)
      DO i2 = 1,ncomp
        fch(i,i2,jx,jy,jz) = fch(i,i2,jx,jy,jz) + mutemp*musurf(ns,i2)*spex_conc
      END DO
      DO is2 = 1,nsurf
        fch(i,is2+ncomp+nexchange,jx,jy,jz) = fch(i,is2+ncomp+nexchange,jx,jy,jz)  &
         + mutemp*musurf(ns,is2+ncomp)*spex_conc
      END DO
    END IF
  END DO
END DO

RETURN
END SUBROUTINE jac_exchange
!***********************************************************
