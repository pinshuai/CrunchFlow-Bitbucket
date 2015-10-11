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
REAL(DP)                                                           :: MeanVtransverse
REAL(DP)                                                           :: vx_T
REAL(DP)                                                           :: qbar_T

INTEGER(I4B)                                                       :: jx
INTEGER(I4B)                                                       :: jy
INTEGER(I4B)                                                       :: jz

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
     
        
      vx = qx(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )
      vy = qy(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )
      vz = qz(jx,jy,jz)/( 0.5d0*(por(jx,jy,jz)+por(jx,jy,jz+1)) )

      qbar = DSQRT( vx*vx + vy*vy + vz*vz )

      IF (qbar /= 0.0) THEN
        dspx(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*qbar +  &
            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )*vx*vx/qbar
          
!!!        dspx(jx,jy,jz) = alft*qbar +  &
!!!            (alfl-alft)*vx*vx/qbar
        
!!!        dspx(jx,jy,jz) = alfl*vx*vx/qbar + alft*vy*vy/qbar
        
      ELSE
        dspx(jx,jy,jz) = 0.0
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        
      vx = qx(jx,jy,jz) /( 0.5d0*(por(jx,jy,jz)+por(jx+1,jy,jz)) )

      
!!      vx_T = qx(jx,jy+1,jz)/( 0.5d0*(por(jx,jy+1,jz)+por(jx+1,jy+1,jz)) )
      
!!      MeanVtransverse = 0.5*(vx+vx_T) 
      MeanVtransverse = 0.0d0
      
      qbar_T = DSQRT( MeanVtransverse*MeanVtransverse)

      

      
      qbar = DSQRT( vx*vx + vy*vy + vz*vz )
      
      IF (qbar /= 0.0) THEN
!!        dspy(jx,jy,jz) =  alft*qbar +  &
!!            (alfl-alft)*qy(jx,jy,jz)*qy(jx,jy,jz)/qbar
          
!!            dspy(jx,jy,jz) = alft*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*qbar +  &
!!            (alfl-alft)*( 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) )*vy*vy/qbar
          
!!             dspy(jx,jy,jz) = alft*qbar +  &
!!            (alfl-alft)*vy*vy/qbar
             dspy(jx,jy,jz) = 0.5d0*(por(jx,jy,jz)+por(jx,jy+1,jz)) * (alfT*qbar_T )
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
