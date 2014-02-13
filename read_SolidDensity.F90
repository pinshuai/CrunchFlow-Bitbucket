!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:03:46
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE read_SolidDensity(nout,nchem)
USE crunchtype
USE CrunchFunctions
USE params
USE strings
USE concentration

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem


!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lsave

!  This routine reads the solid density, or instructions to calculate this from mineral volume fractions

SolidDensityFrom(nchem) = 3               !! Density from mineral volume fractions

REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  lsave = lzs
  IF (ssch == 'soliddensity' .OR. ssch == 'solid_density') THEN            !! SolidDensity keyword found
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      lsave = lzs
      IF (res == 'n') THEN
        SolidDensity(nchem) = DNUM(ssch)
        SolidDensityFrom(nchem) = 1
      ELSE
        IF (ssch == '-ss_ratio' .OR. ssch == 'solidsolutionratio' .OR. ssch == 'solidsolution_ratio' .OR. ssch == 'solid_solution__ratio') THEN
          id = ids + ls
          CALL sschaine(zone,id,iff,ssch,ids,ls)
          IF (ls /= 0) THEN
            lzs=ls
            CALL convan(ssch,lzs,res)
            lsave = lzs
            IF (res == 'n') THEN
              SolidSolutionRatio(nchem) = DNUM(ssch)
              SolidDensityFrom(nchem) = 2
            ELSE
              WRITE(*,*)
              WRITE(*,*) ' A numerical value for the solid:solution ratio'
              WRITE(*,*) '   should follow the parameter flag "SolidSolutionRatio" '
              WRITE(*,*) ssch(1:lsave)
              WRITE(*,5050) nchem
              WRITE(*,*)
              READ(*,*) 
              STOP
            END IF
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' No values given for solid solution ratio'
            WRITE(*,*) ssch(1:lsave)
            WRITE(*,5050) nchem
            WRITE(*,*)
            READ(*,*) 
            STOP
          END IF
        ELSE IF (ssch == 'calculate' .OR. ssch == 'calculatefromminerals' .OR. ssch == 'calculatefrommineral') THEN     !!  Solid density (and solid:solution ratio) to be calculated from mineral volume fractions
          SolidDensityFrom(nchem) = 3
        ELSE
          SolidDensityFrom(nchem) = 3                   !! Default
        END IF
      END IF
    ELSE                                                !! Nothing following the keyword "SolidDensity"
      WRITE(*,*)
      WRITE(*,*) ' Nothing following the keyword "SolidDensity"'
      WRITE(*,5050) nchem
      WRITE(*,*)
      READ(*,*) 
      STOP
    END IF
    GO TO 300                                           !! Exit if the "SolidDensity" keyword has been found
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_SolidDensity
