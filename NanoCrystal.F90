MODULE  NanoCrystal

  USE crunchtype
  USE params

!!  LOGICAL(LGT)                                      ::  
!!  INTEGER(I4B)                                      :: 
!!  REAL(DP), DIMENSION(:), ALLOCATABLE               :: 
!!  CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE  ::

!!  Logicals
  LOGICAL(LGT),DIMENSION(:), ALLOCATABLE              :: CrystalSizeDistribution

!!  Integers
  INTEGER(I4B)                                        :: nCSD
  INTEGER(I4B), DIMENSION(:,:,:,:,:), ALLOCATABLE     :: nCrystal

!!  Reals

  REAL(DP)                                            :: NucleationScaleFactor

  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: sigma
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: radius
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: xtvd
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: aaCSD
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: bbCSD
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: ccCSD
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: uuCSD
  REAL(DP), DIMENSION(:), ALLOCATABLE                 :: rrCSD

  REAL(DP), DIMENSION(:,:), ALLOCATABLE               :: AdjustKeqLog

  REAL(DP), DIMENSION(:,:), ALLOCATABLE               :: AreaCrystalSize

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE             :: rateCrystalSize
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE             :: silogCrystalSize
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE             :: siCrystalSize

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE         :: LinearGrowthRate

!!  Characters


END MODULE NanoCrystal

