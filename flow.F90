MODULE flow
  
USE crunchtype
      
REAL(DP)                                     :: permmaxX
REAL(DP)                                     :: permmaxY
REAL(DP)                                     :: permmaxZ
REAL(DP)                                     :: permminX
REAL(DP)                                     :: permminY
REAL(DP)                                     :: permminZ
REAL(DP)                                     :: maxQx
REAL(DP)                                     :: maxQy
REAL(DP)                                     :: maxQz
REAL(DP)                                     :: x_angle
REAL(DP)                                     :: y_angle
REAL(DP)                                     :: z_angle
REAL(DP)                                     :: SignGravity

LOGICAL(LGT)                                 :: CalculateFlow
LOGICAL(LGT)                                 :: Neumann
LOGICAL(LGT)                                 :: InitializeHydrostatic

INTEGER(I4B)                                 :: infiltration

!fp! block(1);
INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE  :: intbnd
INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE  :: intbndgas
!fp! darray=0;        
!fp! block_xy(1);
REAL(DP), DIMENSION(:,:), ALLOCATABLE        :: qrecharge
!fp! darray=0;

REAL(DP), DIMENSION(:,:,:), ALLOCATABLE      :: gaspump

REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzonex
REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzoney
REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzonez

INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE  :: activecellPressure

INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermx_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermx_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermx_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermx_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermx_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermx_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermy_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermy_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermy_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermy_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermy_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermy_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermz_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermz_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermz_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermz_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermz_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermz_hi

REAL(DP), DIMENSION(:), ALLOCATABLE          :: PressureZone
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: PressureFix
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxPressure_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxPressure_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyPressure_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyPressure_hi
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzPressure_lo
INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzPressure_hi

REAL(DP), DIMENSION(:,:,:),ALLOCATABLE       :: pres
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: harx
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: hary
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: harz
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminx
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminy
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminz
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permx
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permy
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permz
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permxOld
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permyOld
REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permzOld

!!  PETSc arrays for solver

REAL(DP), DIMENSION(:),ALLOCATABLE           :: XvecCrunchP
REAL(DP), DIMENSION(:),ALLOCATABLE           :: BvecCrunchP

END MODULE flow
