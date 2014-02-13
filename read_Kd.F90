SUBROUTINE read_kd(nout,ncomp,ncnt)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
 
IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(OUT)                                   :: ncnt

!  Internal variables and arrays

LOGICAL(LGT)                                                :: speciesnotfound
CHARACTER (LEN=mls), DIMENSION(ncomp)                       :: nam_retard
CHARACTER (LEN=mls)                                         :: namtemp
CHARACTER (LEN=mls)                                          :: dumstring
REAL(DP), DIMENSION(ncomp)                                  :: tempkd

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lspecies
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: lss

speciesnotfound = .true.
REWIND nout

ncnt = 0
100 READ(nout,'(a)',END=300) zone
dumstring = zone(1:13)
CALL majuscules(dumstring,13)
IF (dumstring == 'solid_density') GOTO 100
nlen1 = LEN(zone)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  ncnt = ncnt + 1
  IF (ncnt > ncomp) THEN
    WRITE(*,*)
    WRITE(*,*) ' Too many species for retardation specified'
    WRITE(*,*) ' Number of retarded species should not exceed # of components'
    WRITE(*,*) ' Number of components = ',ncomp
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  lzs=ls
  nam_retard(ncnt) = ssch
  namtemp = nam_retard(ncnt)
  CALL stringlen(namtemp,lspecies)
  
!  Now read the Kd (L/kg)
  
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      tempkd(ncnt) = DNUM(ssch)
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Looking for a number for the Kd in RETARDATION section'
      WRITE(*,*) ' Following primary species ',namtemp(1:lspecies)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No Kd following primary species name in section RETARDATION'
    WRITE(*,*) ' Primary species = ',namtemp(1:lspecies)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
END IF

GO TO 100  ! Keep reading from file until the end of list is found

!  Go here when the list of retarded species is complete--check against primary species list

300 DO kk = 1,ncnt
  namtemp = nam_retard(kk)
  CALL stringlen(namtemp,lss)
  speciesnotfound = .true.
  DO i = 1,ncomp
    IF (ulab(i) == nam_retard(kk)) THEN
      distrib(i) = tempkd(kk)
      speciesnotfound = .false.
    END IF
  END DO
  IF (speciesnotfound) THEN
    WRITE(*,*)
    WRITE(*,*) ' ERROR in RETARDATION BLOCK'
    WRITE(*,*) ' Species associated with distribution coefficient not found in'
    WRITE(*,*) '   primary species list'
    WRITE(*,*) ' Looking for species: ',namtemp(1:lss)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO


END SUBROUTINE read_kd
