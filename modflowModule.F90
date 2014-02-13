MODULE modflowModule
  
  USE crunchtype
  USE params

  INTEGER(I4B)                                                     :: nwells
  INTEGER(I4B)                                                     :: ncnh
  INTEGER(I4B)                                                     :: nrivers
  INTEGER(I4B)                                                     :: ndrains
  INTEGER(I4B)                                                     :: nstress
  INTEGER(I4B)                                                     :: ntModFlow
  INTEGER(I4B)                                                     :: NumModFlowSteps

  REAL(DP)                                                         :: ModFlowCnv
  REAL(DP)                                                         :: timeModFlow

  LOGICAL(LGT)                                                     :: ModFlowRead
  
  REAL(DP), DIMENSION(:), ALLOCATABLE                              :: q
  REAL(DP), DIMENSION(:), ALLOCATABLE                              :: qcnh
!!  REAL(DP), DIMENSION(:,:), ALLOCATABLE                            :: qrecharge
  REAL(DP), DIMENSION(:,:), ALLOCATABLE                            :: qevt
  REAL(DP), DIMENSION(:), ALLOCATABLE                              :: qriver
  REAL(DP), DIMENSION(:), ALLOCATABLE                              :: qdrain

  REAL(DP), DIMENSION(:), ALLOCATABLE                              :: dtModFlow


  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jxWellLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jyWellLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jzWellLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jxDrainLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jyDrainLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jzDrainLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jxRiverLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jyRiverLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jzRiverLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jxHeadLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jyHeadLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: jzHeadLoc
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: WellCondition
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: RiverCondition
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                          :: HeadCondition
  INTEGER(I4B)                                                     :: RechargeCondition

!fp! block_xy(1);
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                        :: lrecharge
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                        :: levt
!fp! darray=0;

  CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE                    :: wtype
  CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE                    :: htype
  CHARACTER (LEN=20), DIMENSION(:), ALLOCATABLE                    :: rtype

!!  LOGICAL(LGT)                                                :: 

!!  REAL(DP)                                                    :: 

!!  CHARACTER (LEN=mls)                                         :: 


END MODULE modflowModule

