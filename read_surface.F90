SUBROUTINE read_surface(nout,ncomp,nkin,nsurf)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nkin
INTEGER(I4B), INTENT(OUT)                                   :: nsurf

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: namtmp
CHARACTER (LEN=mls)                                         :: namsurfsave

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: lsurf
INTEGER(I4B)                                                :: lmin
INTEGER(I4B)                                                :: k

REWIND nout

nsurf = 0

100 READ(nout,'(a)',END=111) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  IF (ssch(1:1) == '>') THEN
    lsurf = ls
    nsurf = nsurf + 1
    IF (nsurf > msurf) THEN
      WRITE(*,*)
      WRITE(*,*) ' Number of surface sites dimensioned too small'
      WRITE(*,*) ' Msurf = ',msurf
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    namsurf(nsurf) = ssch
  ELSE
    GO TO 100
  END IF
  namsurfsave = namsurf(nsurf)

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch /= 'on') THEN
      WRITE(*,*) ' Surface complex should be followed by "on" '
      WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    lzs=ls
    lmin = ls
    CALL stringtype(ssch,lzs,res)
    IF (res /= 'a') THEN
      WRITE(*,*)
      WRITE(*,*) ' Looking for a mineral name, not a number'
      WRITE(*,*) ' In surface complexation block'
      WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    namtmp = ssch
    DO k = 1,nkin
      IF (namtmp == umin(k)) THEN
        ksurf(nsurf) = k
        GO TO 200
      END IF
    END DO
    WRITE(*,*)
    WRITE(*,*) ' Mineral substrate listed in surface complexation block'
    WRITE(*,*) '     not found in minerals list'
    WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
    WRITE(*,*) ' Looking for mineral: ',namtmp(1:lmin)
    WRITE(*,*)
    READ(*,*)
    STOP
    
    200    CONTINUE
    
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Mineral host for surface hydroxyl site must be specified'
    WRITE(*,*) ' In surface complexation block'
    WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch == '-no_edl') THEN
      iedl(nsurf) = 1
    END IF
  ELSE
    GO TO 100
  END IF
  
END IF

GO TO 100

111  CONTINUE

RETURN

END SUBROUTINE read_surface
