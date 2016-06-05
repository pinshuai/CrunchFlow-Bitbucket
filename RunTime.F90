!! CrunchTope 
!! Copyright (c) 2016, Carl Steefel
!! Copyright (c) 2016, The Regents of the University of California, 
!! through Lawrence Berkeley National Laboratory (subject to 
!! receipt of any required approvals from the U.S. Dept. of Energy).  
!! All rights reserved.

!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met: 

!! (1) Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.

!! (2) Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.

!! (3) Neither the name of the University of California, Lawrence
!! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of    
!! its contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.

!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE 

    
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
  LOGICAL(LGT)                                                :: visit
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
  LOGICAL(LGT)                                                :: Duan2006
  LOGICAL(LGT)                                                :: UseBulkMineral
  LOGICAL(LGT)                                                :: Maggi
  LOGICAL(LGT)                                                :: DePaolo
  LOGICAL(LGT)                                                :: SetSurfaceAreaConstant
!!  LOGICAL(LGT)                                                :: UseAqueousMoleFraction
  LOGICAL(LGT)                                                :: ReadGeochemicalConditions
  LOGICAL(LGT)                                                :: ReadGautier
  LOGICAL(LGT)                                                :: Qingyun
  LOGICAL(LGT)                                                :: ForsteriteCapillary
  LOGICAL(LGT)                                                :: CylindricalDivideVolume
  LOGICAL(LGT)                                                :: Benchmark
  LOGICAL(LGT)                                                :: DampRateInLowPorosity

  LOGICAL(LGT)                                                :: Switcheroo
  
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

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE                       :: SkipAdjust

END MODULE runtime
