MODULE isotope

  USE crunchtype
  USE params

  INTEGER(I4B)                                      :: nIsotopePrimary
  INTEGER(I4B)                                      :: nIsotopeMineral
 
!  Permanent arrays, but can be reallocated to smaller size


!! Primary species

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: isotopeRare
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: isotopeCommon

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: IsotopePrimaryRare
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: IsotopePrimaryCommon

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iPointerIsotope

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: nameIsotopeRare
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: nameIsotopeCommon

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: IsotopeReference

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: MoleFractionAqueousRare
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: MoleFractionAqueousCommon

  REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: dMoleFractionAqueousCommon
  REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: dMoleFractionAqueousRare

!! Minerals

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kIsotopeRare
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kIsotopeCommon

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: IsotopeMineralRare
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: IsotopeMineralCommon

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kPointerIsotope

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: nameIsotopeMineralRare
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: nameIsotopeMineralCommon

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: isotopeBackReactionOption

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: MoleFractionMineralRare
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: MoleFractionMineralCommon 

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: PointerToPrimaryIsotope

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: UseAqueousMoleFraction



! Allocatable arrays that are local (not dimensioned over spatial domain)

!  Allocatable arrays dimensioned over spatial domain



END MODULE isotope

