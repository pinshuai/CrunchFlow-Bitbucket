
! Code converted using TO_F90 by Alan Miller
! Date: 2000-08-08  Time: 14:03:41
 
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

SUBROUTINE jacexchange(ncomp,nexchange,nexch_sec,neqn,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE solver
USE transport
USE medium
USE temperature

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                  :: ncomp
INTEGER(I4B), INTENT(IN)                                  :: nexchange
INTEGER(I4B), INTENT(IN)                                  :: nexch_sec
INTEGER(I4B), INTENT(IN)                                  :: neqn
INTEGER(I4B), INTENT(IN)                                  :: jx
INTEGER(I4B), INTENT(IN)                                  :: jy
INTEGER(I4B), INTENT(IN)                                  :: jz

!  Internal variables and arrays

INTEGER(I4B)                                              :: ix
INTEGER(I4B)                                              :: i2
INTEGER(I4B)                                              :: nex
INTEGER(I4B)                                              :: ixcheck

REAL(DP)                                                  :: spex_conc
REAL(DP)                                                  :: mutemp
REAL(DP)                                                  :: exchangetemp
REAL(DP)                                                  :: sum2

fexch = 0.0

IF (iexc == 1 .OR. iexc == 3) THEN
  
  DO nex = 1,nexch_sec
    spex_conc = spex10(nex+nexchange,jx,jy,jz)
    DO ix = 1,nexchange
      IF (muexc(nex,ix+ncomp) /= 0.0) THEN
        mutemp = muexc(nex,ix+ncomp)
        exchangetemp = exchangesites(ix,jx,jy,jz)
        DO i2 = 1,ncomp+nexchange
          fexch(ix+ncomp,i2) = fexch(ix+ncomp,i2) + mutemp*muexc(nex,i2)*spex_conc/  &
            exchangetemp
        END DO
      END IF
    END DO
  END DO

!!        do ix = 1,nexchange
!!          do i2 = 1,ncomp+nexchange
!!            sum2 = 0.0
!!            do nex = 1,nexch_sec
!!              sum2 = sum2 + muexc(nex,ix+ncomp)*muexc(nex,i2)* &            
!!           spex10(nexchange+nex,jx,jy,jz)
!!           end do
!!            fexch(ix+ncomp,i2) = sum2
!!          end do
!!        end do
  
  
ELSE    !                  Vanselow convention

  fweight = 0.0
  
  DO nex = 1,nexch_sec
    spex_conc = aexch(nex)
    ixcheck = ixlink(nex)
    DO ix = 1,nexchange
      IF (muexc(nex,ix+ncomp) /= 0.0) THEN
        mutemp = muexc(nex,ix+ncomp)
        DO i2 = 1,ncomp+nexchange
          fweight(ix+ncomp,i2) = fweight(ix+ncomp,i2) + mutemp*muexc(nex,i2)*spex_conc
        END DO
      END IF
      IF (ixcheck == ix) THEN
        DO i2 = 1,ncomp+nexchange
            fexch(ix+ncomp,i2) = fexch(ix+ncomp,i2) + muexc(nex,i2)*spex_conc
        END DO
      END IF
    END DO
  END DO
  
  
END IF

RETURN
END SUBROUTINE jacexchange
!***********************************************************
