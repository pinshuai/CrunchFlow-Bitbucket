!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:40
 
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

SUBROUTINE read_decay(nout,ncomp,ndecay)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(OUT)                                   :: ndecay

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: namdum

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: idecaytmp
INTEGER(I4B)                                                :: lsave
INTEGER(I4B)                                                :: isotope


ndecay = 0

100 READ(nout,'(a)',END=300) zone

id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  WRITE(*,*) ssch
  DO i = 1,ncomp
    IF (ssch == ulab(i)) THEN
      idecaytmp = i
      GO TO 110
    END IF
  END DO
  WRITE(*,*)
  WRITE(*,*) ' Species undergoing radioactive decay not found in primary species list'
  WRITE(*,*) ' Looking for species: ',ulab(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
  
  110   ndecay = ndecay + 1
  idecay(ndecay) = idecaytmp     ! Pointer to primary species number
  lsave = ls
  isotope = 0
  
!  Now search for the isotopes
  
  120   id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  
  IF (ls /= 0) THEN
    
    isotope = isotope + 1
    decay_label(isotope,ndecay) = ssch
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs = ls
      CALL stringtype(ssch,lzs,res)
      IF (res == 'a') THEN
        namdum = ulab(idecaytmp)
        WRITE(*,*)
        WRITE(*,*) ' Isotope label should be followed by a half-life'
        WRITE(*,*) ' Isotope of   : ',namdum(1:lsave)
        WRITE(*,*) ' Found string : ',ssch(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      ELSE
        half_life(isotope,ndecay) = DNUM(ssch)
      END IF
    ELSE
      namdum = ulab(idecaytmp)
      WRITE(*,*)
      WRITE(*,*) ' Isotope label should be followed by a half-life'
      WRITE(*,*) ' Isotope of: ',namdum(1:lsave)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    GO TO 120   !  Search for another isotope
    
  END IF
  
  nisotope(ndecay) = isotope
  
  GO TO 100     !  Search for another primary species subject to radioactive decay
  
ELSE
  GO TO 100
END IF

300 RETURN
END SUBROUTINE read_decay
