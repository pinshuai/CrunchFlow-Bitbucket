MODULE solver

  USE crunchtype
  USE params

  INTEGER(I4B)                                     :: level
  INTEGER(I4B)                                     :: Gimrtlevel
  REAL(DP)                                         :: GimrtRTOLKSP
  
  CHARACTER (LEN=mls)                              :: GIMRT_SolverMethod
  CHARACTER (LEN=mls)                              :: SolverMethod
  CHARACTER (LEN=mls)                              :: PCMethod
  CHARACTER (LEN=mls)                              :: GIMRT_PCMethod

  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockfl
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockl
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockm
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockr
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockfr
  
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fsurf_local
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fch_local
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fjac_loc
  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE        :: alf
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: aaa
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fexch
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fsurftmp
  REAL(DP), DIMENSION(:),   ALLOCATABLE            :: jac_sat
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: jac_pre
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: jac_preKin
  REAL(DP), DIMENSION(:),   ALLOCATABLE            :: jac_rateFactor

  INTEGER(I4B), DIMENSION(:),    ALLOCATABLE       :: indd

  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: fxmax
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: bb
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: xn
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: fxx

!  Arrays for WATSOLV

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: list
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: lrow
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: levptr
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: row

!  ALLOCATABLE arrays DIMENSIONed over spatial domain

!fp! block_gb(3,1);
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac_D
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac_chg
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fgas
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fsurf
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fch
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjpotncomp
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjpotnsurf
!fp! darray=0;

! Allocatable arrays that are local

  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fjpotncomp_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fjpotnsurf_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fgas_local

! One DIMENSIONal arrays for block tridiagonal solver

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: aah
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: bbh
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: cch
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: yh

  INTEGER(I4B), DIMENSION(:,:),  ALLOCATABLE       :: indexx

! Arrays for PETSc

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: XvecCrunchD
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: BvecCrunchD

 
END MODULE solver
     
