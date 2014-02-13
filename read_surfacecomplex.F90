SUBROUTINE read_surfacecomplex(nout,ks,isolution,  &
    ncomp,nkin,ngas,complexfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ks
INTEGER(I4B), INTENT(IN)                                    :: isolution
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nkin
INTEGER(I4B), INTENT(IN)                                    :: ngas
LOGICAL(LGT), INTENT(IN OUT)                                :: complexfound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: namtemp

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: lsurf

CALL stringlen(namtemp,lsurf)

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  IF (ssch == namsurf(ks)) THEN
    complexfound = .true.     ! Surface complex found
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
!  Read the concentration
    site_density(ks,isolution) = DNUM(ssch)
  ELSE
!  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A number should follow surface complex name'
    WRITE(*,*) ' For surface complex ',namtemp(1:lsurf)
    WRITE(*,5050) isolution
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No site density or specific surface area given '
  WRITE(*,*) ' For surface complex ',namtemp(1:lsurf)
  WRITE(*,5050) isolution
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!id = ids + ls
!CALL sschaine(zone,id,iff,ssch,ids,ls)

!IF(ls /= 0) THEN
!  lzs=ls
!  CALL convan(ssch,lzs,res)
!  IF (res == 'n') THEN
!    specific(ks,isolution) = DNUM(ssch)
!  ELSE
!  An ascii string--so bag it.
!    WRITE(*,*)
!    WRITE(*,*) ' A number should follow surface complex name'
!    WRITE(*,*) ' For surface complex ',namtemp(1:lsurf)
!    WRITE(*,5050) isolution
!    WRITE(*,*)
!    STOP
!  END IF
  
!ELSE      !  No trailing string found
!  WRITE(*,*)
!  WRITE(*,*) ' No specific surface area given '
!  WRITE(*,*) ' For surface complex ',namtemp(1:lsurf)
!  WRITE(*,5050) isolution
!  WRITE(*,*)
!  STOP
!END IF

5050 FORMAT(1X,'Geochemical condition number: ',i2,' found')

300  RETURN
END SUBROUTINE read_surfacecomplex
