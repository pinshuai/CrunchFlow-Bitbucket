!************** (C) COPYRIGHT 1995 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:01
 
!                      All Rights Reserved

!  OSRT (Operator Splitting Reactive Transport) IS PROVIDED "AS IS"
!  AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED. THE USER ASSUMES ALL RISKS
!  OF USING OSRT. THERE IS NO CLAIM OF THE MERCHANTABILITY OR FITNESS
!  FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

SUBROUTINE timestep(nx,ny,nz,delt,dtold,ttol,tstep,dtmax,ikmast)
USE crunchtype
USE params
USE concentration

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz
INTEGER(I4B), INTENT(IN)                                        :: ikmast

REAL(DP), INTENT(IN OUT)                                        :: delt
REAL(DP), INTENT(IN OUT)                                        :: dtold
REAL(DP), INTENT(IN)                                            :: ttol
REAL(DP), INTENT(IN)                                            :: tstep
REAL(DP), INTENT(IN)                                            :: dtmax

!  Internal variables and arrays

REAL(DP), PARAMETER                                             :: theta=0.8
REAL(DP)                                                        :: deltmax
REAL(DP)                                                        :: stepmax
REAL(DP)                                                        :: dudtp
REAL(DP)                                                        :: dudt
REAL(DP)                                                        :: rnum
REAL(DP)                                                        :: rden
REAL(DP)                                                        :: dtemp
REAL(DP)                                                        :: maxtstep

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz

maxtstep = MIN(dtmax,tstep)
deltmax = 2.0*delt
stepmax = MIN(deltmax,maxtstep)

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      dudtp = (sp(ikmast,jx,jy,jz) - spno2(jx,jy,jz))/delt  
!! Monitor "master" variable, generally pH or redox species
      dudt = (spno2(jx,jy,jz) - spnno2(jx,jy,jz))/dtold        
!!  spno2 is concentration from last time step, spnno2 from the time step before that
      rnum = 2.0D0*theta*ttol*(spno2(jx,jy,jz)+sp(ikmast,jx,jy,jz))  !! ttol is time tolerance
      rden = (2.0D0/(delt+dtold))*(dudtp-dudt)  
!! "rden" gives the second derivative of master concentration change with respect to time (d2C/dt2)
      IF (ABS(rden) < 1.0E-30) GO TO 20
        dtemp = DSQRT(DABS(rnum/rden))
        IF(dtemp/delt > 2.0D0) THEN
          dtemp = 2.0D0*delt  !!  Only increase by a factor of 2 at most
        END IF
        IF (dtemp > stepmax) GO TO 15   !! If greater than "stepmax", use "stepmax"
        stepmax = dtemp
        15       CONTINUE
        20       CONTINUE
     END DO
  END DO
END DO
dtold = delt
delt = stepmax

RETURN
END SUBROUTINE timestep
!  *******************************************************
