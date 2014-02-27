!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:05:39
 
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

SUBROUTINE read_mineral(nout,k,isolution,ncomp,nspec,nrct,ngas,mineralfound)
USE crunchtype
USE params
USE mineral
USE concentration
USE strings
USE CrunchFunctions, ONLY: dnum

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: k
INTEGER(I4B), INTENT(IN)                                    :: isolution
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nrct
INTEGER(I4B), INTENT(IN)                                    :: ngas
LOGICAL(LGT), INTENT(IN OUT)                                :: mineralfound

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lmin
INTEGER(I4B)                                                :: lcond

CHARACTER (LEN=mls)                                         :: minlabel
CHARACTER (LEN=mls)                                         :: labeltemp

minlabel = umin(k)
CALL stringlen(minlabel,lmin)


labeltemp = condlabel(isolution)
CALL stringlen(labeltemp,lcond)

REWIND nout

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
!      call majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
  IF (res /= 'a') THEN
    WRITE(*,*)
    WRITE(*,*) ' Geochemical input should start with a string'
    WRITE(*,*) ssch(1:ls)
    WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF
IF (ssch == umin(k)) THEN
  mineralfound = .true.     ! mineral found
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
    volin(k,isolution) = DNUM(ssch)
  ELSE
    IF (ssch=='mole' .OR. ssch=='moles') THEN

      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
!       Read the concentration
          MineralMoles(k,isolution) = DNUM(ssch)
        ELSE
!         An ascii string--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' A number should follow the mineral name'
          WRITE(*,*) ' For mineral ', minlabel(1:lmin)
          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
          READ(*,*)
          STOP
        END IF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No mineral concentration given '
        WRITE(*,*) ' For mineral ', minlabel(1:lmin)
        WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
!     An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A number should follow the mineral name'
      WRITE(*,*) ' For mineral ', minlabel(1:lmin)
      WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
      READ(*,*)
      STOP
    END IF
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*) ' No volume fraction given '
  WRITE(*,*) ' For mineral ', minlabel(1:lmin)
  WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
  WRITE(*,*) ' Assuming mineral initially not present'
  WRITE(*,*)
  volin(k,isolution) = 0.0
END IF
IF (volin(k,isolution) >= 1.0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Volume fraction of mineral > 1.0'
  WRITE(*,*) ' For mineral ', minlabel(1:lmin)
  WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)
  WRITE(*,*)
!!  READ(*,*)
!!  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)

IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
    areain(k,isolution) = DNUM(ssch)
    iarea(k,isolution) = 0 
  ELSE
!  An ascii string
    IF (ssch == 'bulk_surface_area' .OR. ssch == 'bsa') THEN
      iarea(k,isolution) = 0 
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF (ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          areain(k,isolution) = DNUM(ssch) 
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' No number following specification of bulk_surface_area'
          WRITE(*,*) ' For mineral ', minlabel(1:lmin)
          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond) 
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No number following specification of bulk_surface_area'
        WRITE(*,*) ' For mineral ', minlabel(1:lmin)
        WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)  
        WRITE(*,*)
        READ(*,*)
        STOP
      ENDIF   
    ELSE IF (ssch == 'specific_surface_area' .OR. ssch == 'ssa') THEN
      iarea(k,isolution) = 1 
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF (ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          specific(k,isolution) = DNUM(ssch)  
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' No number following specification of specific_surface_area'
          WRITE(*,*) ' For mineral ', minlabel(1:lmin)
          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)   
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No number following specification of specific_surface_area'
        WRITE(*,*) ' For mineral ', minlabel(1:lmin)
        WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)   
        WRITE(*,*)
        READ(*,*)
        STOP
      ENDIF 
!!  Now look for a "threshold" volume fraction that is used when the initial volume fraction = 0.  This value will be used to calculate a bulk surface area given a specific
!!  surface area, until the actual computed volume fraction of the secondary mineral reaches this value.  At this stage, the actual computed volume fraction is used.
!!  Do this ONLY if volin(k,isolution) = 0
      IF (volin(k,isolution) == 0.0d0 .AND. MineralMoles(k,isolution)==0.0d0) THEN
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF (ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            voltemp(k,isolution) = DNUM(ssch)  
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Specification of a specific surface area should be followed by a threshold volume fraction where initial volume fraction = 0'
            WRITE(*,*) ' For mineral ', minlabel(1:lmin)
            WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)   
            WRITE(*,*)
            READ(*,*)
            STOP
          ENDIF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Specification of a specific surface area should be followed by a threshold volume fraction where initial volume fraction = 0'
          WRITE(*,*) ' For mineral ', minlabel(1:lmin)
          WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)   
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIF 
      END IF

    END IF  
  END IF
  
ELSE      !  No trailing string found
!!  WRITE(*,*)
!!  WRITE(*,*) '    No mineral surface area given '
!!        WRITE(*,*) ' For mineral ', minlabel(1:lmin)
!!        WRITE(*,*) ' In geochemical condition: ', labeltemp(1:lcond)  
!   write(*,*) ' Assuming 100 m**2/m**3'
!!  WRITE(*,*)
!!  READ(*,*)
!!  STOP
   areain(k,isolution) = 1.0
END IF

5050 FORMAT(1X,'Geochemical condition number ',i2,' found')

300  RETURN
END SUBROUTINE read_mineral
