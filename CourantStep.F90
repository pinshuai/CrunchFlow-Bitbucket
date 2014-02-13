!  ************ Subroutine TIMESTEP ***********
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-09-18  Time: 10:34:23
 
!  Purpose is to check courant condition (courfactor) for all MODFLOW flows and sources
!  Written by Steve Yabusaki and Ashok Chilakapati 6/2000
!  *****************************************************************************

SUBROUTINE CourantStep(nx,ny,nz,dtmaxcour)

USE crunchtype
USE runtime
USE medium
USE transport
USE flow
USE modflowModule

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

REAL(DP), INTENT(OUT)                                    :: dtmaxcour

!  Internal variables and arrays

REAL(DP)                                                 :: dtmin,dt,cnv

INTEGER(I4B)                                             :: jx,jy,jz
INTEGER(I4B)                                             :: i

cnv = ModFlowCnv  !  Converts from MODFLOW time units to years (the CRUNCH time unit)

!   calculations for maximum time step based on courant condition

dtmin=1.e30
DO jz=1,nz
  DO jy=1,ny
    DO jx=1,nx
!fp! if_onproc( {#expr# dzz(jx,jy,jz) #});
      IF(qx(jx,jy,jz) == 0.0)THEN
        dt= 1.e30
      ELSE
        dt= ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*por(jx,jy,jz)/qx(jx,jy,jz))/cnv
!fp! end_onproc();
      END IF
      dtmin=MIN(dt,dtmin)
      IF(qy(jx,jy,jz) == 0.0)THEN
        dt=1.e30
      ELSE
        dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*por(jx,jy,jz)/qy(jx,jy,jz))/cnv
      END IF
      dtmin=MIN(dt,dtmin)
    END DO
  END DO
END DO
WRITE(*,*)'qxy:',dtmin

!   wells

DO i = 1, nwells
  IF(q(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jxWellLoc(i),jyWellLoc(i),jzWellLoc(i)) #});
    jx = jxWellLoc(i)
    jy = jyWellLoc(i)
    jz = jzWellLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/q(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'well:',dtmin

!   constant head cells
DO i = 1, ncnh
  IF(qcnh(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jxHeadLoc(i),jyHeadLoc(i),jzHeadLoc(i)) #});
    jx = jxHeadLoc(i)
    jy = jyHeadLoc(i)
    jz = jzHeadLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/qcnh(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'cnh:',dtmin

!   recharge-evapotranspiration
DO jy=1,ny
  DO jx=1,nx
    IF(qrecharge(jx,jy)+qevt(jx,jy) == 0.0)THEN
      dt=1.e30
    ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
      dt=ABS( courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,1)*por(jx,jy,1) /  &
          (qrecharge(jx,jy) + qevt(jx,jy)) )/cnv
!fp! end_onproc();
    END IF
    dtmin=MIN(dt,dtmin)
  END DO
END DO

!   river
DO i=1,nrivers
  IF(qriver(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
    jx = jxRiverLoc(i)
    jy = jyRiverLoc(i)
    jz = jzRiverLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/qriver(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'river:',dtmin

dtmaxcour = dtmin

RETURN
END SUBROUTINE CourantStep
