MODULE temperature

  USE crunchtype
  USE params

  CHARACTER (LEN=mls)                              :: Tfile

  INTEGER(I4B)                                     :: jtemp
  INTEGER(I4B)                                     :: ntemp
  INTEGER(I4B)                                     :: TPointer
  
  LOGICAL(LGT)                                     :: RunIsothermal
    
  REAL(DP)                                         :: tinit
  REAL(DP)                                         :: tgrad

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: tempcond
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: rocond

!  Allocatable arrays dimensioned over spatial domain

!fp! block(1);
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: t
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: told
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ro
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: roOld

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: rogas
!fp! darray=0;

END MODULE temperature
