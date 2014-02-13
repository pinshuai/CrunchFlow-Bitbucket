SUBROUTINE dispersivity(nx,ny,nz)
USE crunchtype
USE transport
USE medium, ONLY: por
USE concentration, ONLY: jinit

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                           :: nx
INTEGER(I4B), INTENT(IN)                                           :: ny
INTEGER(I4B), INTENT(IN)                                           :: nz

!  Internal variables and arrays

REAL(DP)                                                           :: qbar
REAL(DP)                                                           :: vx
REAL(DP)                                                           :: vy
REAL(DP)                                                           :: vz

INTEGER(I4B)                                                       :: jx
INTEGER(I4B)                                                       :: jy
INTEGER(I4B)                                                       :: jz


DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx-1

!!     IF (jinit(jx,jy,jz) == 1) THEN
!!       alfL = 0.01
!!       alfT = 0.001
!!     ELSE IF (jinit(jx,jy,jz) == 2) THEN
!!       alfL = 0.06
!!       alfT = 0.006
!!     ELSE
!!       write(*,*) ' Should not be any other conditions here'
!!       write(*,*)
!!       stop
!!     END IF

!!     IF (jinit(jx,jy,jz) == 1 ) THEN
!!       alfL = 10.0d0
!!       alfT = 1.0d0
!!     ELSE IF (jinit(jx,jy,jz) == 2) THEN
!!       alfL = 60.0d0
!!       alfT = 6.0d0
!!     ELSE
!!       alfL = 0.5d0*(60.0d0+10.0d0)
!!       alfT = 0.5d0*(6.0d0+1.0d0)
!!     END IF

!!       write(*,*) jx,alfL
 
!!     IF (jx > 100 .and. jx < 200) then
!!     IF (jx > 400 .and. jx < 440) then
!!          alfL = 60.0
!!          alfT = 6
!!      ELSE
!!          alfL = 10.0
!!          alfT = 1
!!      ENDIF

!!      IF (jx == 400 .OR. jx == 440) THEN
!!        alfL = 0.5d0*(60+10)
!!      END IF

!!      IF (jx == 100 .OR. jx == 200) THEN
!!        alfL = 0.5d0*(0.07)
!!      END IF

      vx = qx(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )
      vy = qy(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )
      vz = qz(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx,jy,jz+1)) )

      qbar = DSQRT( vx*vx + vy*vy + vz*vz )

      IF (qbar /= 0.0) THEN
        dspx(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*qbar +  &
            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*vx*vx/qbar
      ELSE
        dspx(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      qbar = SQRT( qx(jx,jy,jz)*qx(jx,jy,jz) + qy(jx,jy,jz)*qy(jx,jy,jz)  &
          + qz(jx,jy,jz)*qz(jx,jy,jz) )
      IF (qbar /= 0.0) THEN
        dspy(jx,jy,jz) =  alft*qbar +  &
            (alfl-alft)*qy(jx,jy,jz)*qy(jx,jy,jz)/qbar
      ELSE
        dspy(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      qbar = SQRT( qx(jx,jy,jz)*qx(jx,jy,jz) + qy(jx,jy,jz)*qy(jx,jy,jz)  &
          + qz(jx,jy,jz)*qz(jx,jy,jz) )
      IF (qbar /= 0.0) THEN
        dspz(jx,jy,jz) = alft*qbar +  &
            (alfl-alft)*qz(jx,jy,jz)*qz(jx,jy,jz)/qbar
      ELSE
        dspz(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE dispersivity
