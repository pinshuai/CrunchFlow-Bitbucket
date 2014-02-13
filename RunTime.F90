MODULE runtime
  
  USE crunchtype
  USE params
  
  INTEGER(I4B)                                                :: iterat
  LOGICAL(LGT)                                                :: UseDissolutionOnly
  LOGICAL(LGT)                                                :: checkSaturation

  LOGICAL(LGT)                                                :: gimrt
  LOGICAL(LGT)                                                :: os3d
  LOGICAL(LGT)                                                :: petscon
  LOGICAL(LGT)                                                :: RunToSteady
  LOGICAL(LGT)                                                :: giambalvo
  LOGICAL(LGT)                                                :: KateMaher
  LOGICAL(LGT)                                                :: Inagaki
  LOGICAL(LGT)                                                :: Inagaki2
  LOGICAL(LGT)                                                :: HellmannRateLaw
  LOGICAL(LGT)                                                :: SilicaRateLaw
  LOGICAL(LGT)                                                :: BurchRateLaw
  LOGICAL(LGT)                                                :: OelkersRateLaw
  LOGICAL(LGT)                                                :: xtool
  LOGICAL(LGT)                                                :: tecplot
  LOGICAL(LGT)                                                :: originlab
  LOGICAL(LGT)                                                :: xmgr
  LOGICAL(LGT)                                                :: kaleidagraph
  LOGICAL(LGT)                                                :: nview
  LOGICAL(LGT)                                                :: AppendRestart
  LOGICAL(LGT)                                                :: Cylindrical
  LOGICAL(LGT)                                                :: Rectangular
  LOGICAL(LGT)                                                :: Spherical
  LOGICAL(LGT)                                                :: modflow
  LOGICAL(LGT)                                                :: os3dpetsc 
  LOGICAL(LGT)                                                :: RunningPest
  LOGICAL(LGT)                                                :: CreatePestInstructionFile
  LOGICAL(LGT)                                                :: MakeMovie
  LOGICAL(LGT)                                                :: HanfordStrontium
  LOGICAL(LGT)                                                :: H2Opresent
  LOGICAL(LGT)                                                :: JennyDruhan
  LOGICAL(LGT)                                                :: JennyFirstOrder
  LOGICAL(LGT)                                                :: Duan
  LOGICAL(LGT)                                                :: UseBulkMineral
  LOGICAL(LGT)                                                :: Maggi
  LOGICAL(LGT)                                                :: DePaolo
  LOGICAL(LGT)                                                :: SetSurfaceAreaConstant
!!  LOGICAL(LGT)                                                :: UseAqueousMoleFraction
  LOGICAL(LGT)                                                :: ReadGeochemicalConditions
  LOGICAL(LGT)                                                :: CylindricalDivideVolume
  LOGICAL(LGT)                                                :: Benchmark
  LOGICAL(LGT)                                                :: DampRateInLowPorosity
  
  REAL(DP)                                                    :: PorosityDamp
  REAL(DP)                                                    :: OutputTimeScale
  REAL(DP)                                                    :: OutputDistanceScale
  REAL(DP)                                                    :: steadytol
  REAL(DP), DIMENSION(:), ALLOCATABLE                         :: prtint
  REAL(DP), DIMENSION(:), ALLOCATABLE                         :: OutputTime
  REAL(DP)                                                    :: voltol
  REAL(DP)                                                    :: LagSurface
  REAL(DP)                                                    :: ResidualTolerance


  REAL(DP)                                                    :: courfactor

  INTEGER(I4B)                                                :: irestart
  INTEGER(I4B)                                                :: ihindmarsh
  INTEGER(I4B)                                                :: nchem
  INTEGER(I4B)                                                :: irecharge
  INTEGER(I4B)                                                :: ScreenInterval
  INTEGER(I4B)                                                :: NumOutputTimes
  INTEGER(I4B)                                                :: OutputTimeCounter
  INTEGER(I4B)                                                :: ikTracer
  INTEGER(I4B)                                                :: ncounter

  CHARACTER (LEN=mls)                                         :: restartfile

  CHARACTER (LEN=mls)                                         :: master
  CHARACTER (LEN=mls)                                         :: OutputTimeUnits
  CHARACTER (LEN=mls)                                         :: MODFLOWfile
  CHARACTER (LEN=mls)                                         :: DensityModule
  CHARACTER (LEN=mls)                                         :: RestartOutputFile
  CHARACTER (LEN=mls)                                         :: PestExchangeOutputFile
  CHARACTER (LEN=mls)                                         :: PestSurfaceOutputFile

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: InputFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: TimeSeriesFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: TimeSeriesUnits
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: TimeSeriesSpecies

INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE                   :: IntegerArray
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                       :: RealArray
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: FlexibleArray 

END MODULE runtime
