MODULE medium

  USE crunchtype

  INTEGER(I4B)                                    :: ierode
  INTEGER(I4B)                                    :: isaturate

  REAL(DP)                                        :: constantpor 
  REAL(DP)                                        :: FixSaturation
  REAL(DP)                                        :: MinimumPorosity

  REAL(DP), DIMENSION(:), ALLOCATABLE             :: porcond
  REAL(DP), DIMENSION(:), ALLOCATABLE             :: SaturationCond
  REAL(DP), DIMENSION(:), ALLOCATABLE             :: AqueousToBulkCond

! Allocatable arrays dimensioned over spatial domain

  REAL(DP), dimension(:,:,:), allocatable         :: porin
  REAL(DP), dimension(:,:,:), allocatable         :: por
  REAL(DP), dimension(:,:,:), allocatable         :: porOld
  REAL(DP), dimension(:), allocatable             :: x
  REAL(DP), dimension(:), allocatable             :: y
  REAL(DP), dimension(:), allocatable             :: z
  REAL(DP), dimension(:), allocatable             :: dxx
  REAL(DP), dimension(:), allocatable             :: dyy
  REAL(DP), dimension(:,:,:), allocatable         :: dzz
  REAL(DP), dimension(:,:,:), allocatable         :: dxy


END MODULE medium
