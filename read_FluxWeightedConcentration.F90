
SUBROUTINE read_FluxWeightedConcentration(nout,ncomp,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
USE runtime

IMPLICIT NONE

INTERFACE
  SUBROUTINE read_multstring(nout,lchar,parchar,parfind, stringarray,lenarray,section)
  USE crunchtype
  USE params
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(OUT)                                   :: lchar
  CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
  CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
  CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: stringarray
  INTEGER(I4B), INTENT(OUT)                                   :: lenarray
  CHARACTER (LEN=mls), INTENT(IN)                             :: section
  END SUBROUTINE read_multstring
END INTERFACE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: lchar
INTEGER(I4B)                                                :: lenarray
INTEGER(I4B)                                                :: ll
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: nlength

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: workchar1

!!INTEGER(I4B)                                               :: lchar
CHARACTER (LEN=mls)                                        :: parchar
CHARACTER (LEN=mls)                                        :: parfind
!!CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE             :: stringarray
!!INTEGER(I4B)                                               :: lenarray
CHARACTER (LEN=mls)                                        :: section
CHARACTER (LEN=mls)                                        :: dumstring

nFluxWeightedConcentrationFile = 0

IF (ALLOCATED(FluxWeightedConcentrationFile)) THEN
  DEALLOCATE(FluxWeightedConcentrationFile)
END IF
ALLOCATE(FluxWeightedConcentrationFile(1))

IF (.NOT. ALLOCATED(stringarray)) THEN
  ALLOCATE(stringarray(100))
END IF

IF (ALLOCATED(jxFluxWeightedConcentration_lo)) THEN
  DEALLOCATE(jxFluxWeightedConcentration_lo)
END IF
ALLOCATE(jxFluxWeightedConcentration_lo(100))
IF (ALLOCATED(jxFluxWeightedConcentration_hi)) THEN
  DEALLOCATE(jxFluxWeightedConcentration_hi)
END IF
ALLOCATE(jxFluxWeightedConcentration_hi(100))
IF (ALLOCATED(jyFluxWeightedConcentration_lo)) THEN
  DEALLOCATE(jyFluxWeightedConcentration_lo)
END IF
ALLOCATE(jyFluxWeightedConcentration_lo(100))
IF (ALLOCATED(jyFluxWeightedConcentration_hi)) THEN
  DEALLOCATE(jyFluxWeightedConcentration_hi)
END IF
ALLOCATE(jyFluxWeightedConcentration_hi(100))
IF (ALLOCATED(jzFluxWeightedConcentration_lo)) THEN
  DEALLOCATE(jzFluxWeightedConcentration_lo)
END IF
ALLOCATE(jzFluxWeightedConcentration_lo(100))
IF (ALLOCATED(jzFluxWeightedConcentration_hi)) THEN
  DEALLOCATE(jzFluxWeightedConcentration_hi)
END IF
ALLOCATE(jzFluxWeightedConcentration_hi(100))


nxyz = nx*ny*nz

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check to see if initial substring is "WriteFluxWeightedConcentration" (could be more than 1)
  
  IF (ssch == 'writefluxweightedconcentration') THEN
    nFluxWeightedConcentrationFile = nFluxWeightedConcentrationFile + 1
  ELSE
    GO TO 10
  END IF
ELSE         ! No string found
  GO TO 10
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  i = size(FluxWeightedConcentrationFile,1)
  ALLOCATE(workchar1(i))
  workchar1 = FluxWeightedConcentrationFile
  DEALLOCATE(FluxWeightedConcentrationFile)
  ALLOCATE(FluxWeightedConcentrationFile(nFluxWeightedConcentrationFile))
  IF(nFluxWeightedConcentrationFile /= 0) FluxWeightedConcentrationFile(1:nFluxWeightedConcentrationFile-1) = workchar1(1:nFluxWeightedConcentrationFile-1)
  DEALLOCATE(workchar1)
  FluxWeightedConcentrationFile(nFluxWeightedConcentrationFile) = ssch
ELSE
  IF (nFluxWeightedConcentrationFile == 1) THEN
    FluxWeightedConcentrationFile(1) = 'FluxWeightedConcentration.out'
  ELSE
    WRITE(*,*) 
    WRITE(*,*) ' File name for FluxWeightedConcentrationFile required when more than one is used'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (nxyz == 1) THEN
  jxFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = 1
  jyFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = 1
  jzFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = 1
  jxFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = 1
  jyFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = 1
  jzFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jxFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jxFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jxFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = jxFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for FluxWeightedConcentrationFile'
  WRITE(*,*) ' Flux weighted concentration zone ',nFluxWeightedConcentrationFile
  WRITE(*,*)
  STOP
END IF

!!!  Got to here

IF (ny == 1) THEN
  jyFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = 1
  jzFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = 1
  jyFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = 1
  jzFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jyFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jyFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jyFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = jyFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for FluxWeightedConcentrationFile'
  WRITE(*,*) ' FluxWeightedConcentration zone ',nFluxWeightedConcentrationFile
  WRITE(*,*)
  STOP
END IF

!  Now it should be the Z coordinate

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jzFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jzFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of FluxWeightedConcentrationFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jzFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile) = jzFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No Z grid location given for FluxWeightedConcentrationFile'
  WRITE(*,*) ' FluxWeightedConcentration zone ',nFluxWeightedConcentrationFile
  WRITE(*,*)
  STOP
END IF

!! Go to next line to see if there is another FluxWeightedConcentration file
500 GO TO 10

1000  CONTINUE

REWIND(NOUT)

!! Now read the species names to be written to the FluxWeightedConcentration file

parchar = 'FluxWeightedConcentrationSpecies'
parfind = ' '
lenarray = 0
nplotFluxWeightedConcentration = 0
section = 'OUTPUT'
  
CALL read_multstring(nout,lchar,parchar,parfind,stringarray,lenarray,section)
IF (parfind == 'FluxWeightedConcentrationSpecies' .OR. parfind == 'fluxweightedconcentrationspecies') THEN   

!!  Check to see that strings match species names
      
  nplotFluxWeightedConcentration = lenarray

  IF (ALLOCATED(iplotFluxWeightedConcentration)) THEN
    DEALLOCATE(iplotFluxWeightedConcentration)
  END IF
  ALLOCATE(iplotFluxWeightedConcentration(nplotFluxWeightedConcentration))

  DO ll = 1,nplotFluxWeightedConcentration
    iplotFluxWeightedConcentration(ll) = 0
    IF (stringarray(ll) == 'pH' .OR. stringarray(ll) == 'ph') THEN
      WRITE(*,*) '       ERROR '
      WRITE(*,*) ' Fluxes should be written in terms of total concentrations (not pH)'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    DO ik = 1,ncomp
      IF (ulab(ik) == stringarray(ll)) THEN
        iplotFluxWeightedConcentration(ll) = ik
      END IF
    END DO
    IF (iplotFluxWeightedConcentration(ll) == 0) THEN
      dumstring = stringarray(ll)
      CALL stringlen(dumstring,nlength)
      WRITE(*,*)
      WRITE(*,*) ' Flux weighted concentration species not found in list: ',dumstring(1:nlength)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END DO

ELSE

  nplotFluxWeightedConcentration = 0

END IF
  
IF (ALLOCATED(FluxWeightedConcentrationSpecies)) THEN
  DEALLOCATE(FluxWeightedConcentrationSpecies)
END IF
ALLOCATE(FluxWeightedConcentrationSpecies(nplotFluxWeightedConcentration))
FluxWeightedConcentrationSpecies = stringarray(1:nplotFluxWeightedConcentration)

!!DEALLOCATE(stringarray)

RETURN
END SUBROUTINE read_FluxWeightedConcentration
