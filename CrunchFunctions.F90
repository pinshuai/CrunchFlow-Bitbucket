MODULE CrunchFunctions

CONTAINS

!!  *******************************************

!!  PURE FUNCTION jnum (string)
!!  USE CrunchType
!!  USE params
!!  IMPLICIT NONE
!!  CHARACTER (LEN=mls), INTENT(IN)          :: string
!!  INTEGER(I4B)                           :: jnum
!!  READ(string, *) jnum
!!  END FUNCTION jnum
!! **********************************************
!!  PURE FUNCTION dnum (string)
!!  USE CrunchType
!!  USE params
!!  IMPLICIT NONE
!!  CHARACTER (LEN=mls), INTENT(IN)          :: string
!!  REAL(DP)                               :: dnum
!!  READ(string, *) dnum
!!  END FUNCTION dnum
!! **********************************************
  PURE FUNCTION IntegerToCharacter (i)
  USE CrunchType
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)               :: i
  CHARACTER (LEN=mls)                    :: IntegerToCharacter
  CHARACTER (LEN=mls)                    :: dumstring
  WRITE(dumstring,*) i
  IntegerToCharacter = ADJUSTL(dumstring)
  END FUNCTION IntegerToCharacter

!! **********************************************
  PURE FUNCTION RealToCharacter (x)
  USE CrunchType
  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN)                   :: X
  CHARACTER (LEN=mls)                    :: RealToCharacter
  CHARACTER (LEN=mls)                    :: dumstring
  WRITE(dumstring,*) x
  RealToCharacter = ADJUSTL(dumstring)
  END FUNCTION RealToCharacter
!! **********************************************
PURE FUNCTION ArithmeticMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: ArithmeticMean
ArithmeticMean = 0.5d0*(value1+value2)
END FUNCTION ArithmeticMean
!! **********************************************
PURE FUNCTION HarmonicMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: HarmonicMean
REAL(DP)                                      :: denominator
denominator = value1 + value2
IF (denominator /= 0.0d0) THEN
  HarmonicMean = 2.0d0*value1*value2/denominator
ELSE
  HarmonicMean = 0.0d0
END IF
END FUNCTION HarmonicMean
!! **********************************************
PURE FUNCTION GeometricMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: GeometricMean
GeometricMean = DSQRT(value1*value2)
END FUNCTION GeometricMean
!! **********************************************
PURE FUNCTION GetPrimarySpeciesNumber (ncomp,dumstring)
USE CrunchType
USE params, ONLY: mls
USE concentration,ONLY: ulab
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)               :: ncomp
CHARACTER (LEN=mls), INTENT(IN)        :: dumstring
INTEGER(I4B)                           :: GetPrimarySpeciesNumber
INTEGER(I4B)                           :: i
INTEGER(I4B)                           :: ls
GetPrimarySpeciesNumber = 0
DO i = 1,ncomp
  IF (dumstring == ulab(i)) THEN
    GetPrimarySpeciesNumber = i
  END IF
END DO
END FUNCTION GetPrimarySpeciesNumber

PURE FUNCTION imaxloc(arr)
USE crunchtype
REAL(DP), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B) :: imaxloc
INTEGER(I4B), DIMENSION(1) :: imax

imax=maxloc(arr(:))
imaxloc=imax(1)
END FUNCTION imaxloc

PURE FUNCTION outerprod(a,b)
USE crunchtype
REAL(DP), DIMENSION(:), intent(in)                         :: a,b
REAL(DP), DIMENSION(size(a),size(b))                       :: outerprod

outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod

END MODULE CrunchFunctions