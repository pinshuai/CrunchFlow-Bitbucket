MODULE ReadFlow

  USE crunchtype
  USE params

  CHARACTER (LEN=mls)                                        :: NuftFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE             :: NuftVariable

  LOGICAL(LGT)                                               :: NuftRead
  LOGICAL(LGT)                                               :: ReadNuft
  LOGICAL(LGT)                                               :: AlreadyUnsaturated

  REAL(DP)                                                   :: timeNuft
  REAL(DP)                                                   :: timeNuftOld

  INTEGER(I4B)                                               :: nxNuft
  INTEGER(I4B)                                               :: nyNuft
  INTEGER(I4B)                                               :: nzNuft
  INTEGER(I4B)                                               :: iunitNuft
  INTEGER(I4B)                                               :: NumNuftVariables
  INTEGER(I4B)                                               :: NumSourceTerms
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jxNuftSource
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jyNuftSource
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jzNuftSource
  INTEGER(I4B)                                               :: NumNuftSteps
  INTEGER(I4B)                                               :: TotNuftSteps
  REAL(DP)                                                   :: NuftWeightNew
  REAL(DP)                                                   :: NuftWeightOld

  LOGICAL(LGT)                                               :: ReadNuftLiqFlux
  LOGICAL(LGT)                                               :: ReadNuftSaturation
  LOGICAL(LGT)                                               :: ReadNuftPorosity
  LOGICAL(LGT)                                               :: ReadNuftTemperature
  LOGICAL(LGT)                                               :: ReadNuftLiqDensity
  LOGICAL(LGT)                                               :: ReadNuftGasFlux
  LOGICAL(LGT)                                               :: ReadNuftGasDensity

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qyNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qyNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxgasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qygasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzgasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxgasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qygasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzgasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: satliqNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: satliqNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qgNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qgNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: roNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: roNuftOld


END MODULE ReadFlow