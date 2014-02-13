!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:03:17
 
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

SUBROUTINE read_exchange(nout,ncomp,nexchange,data1,nexch_sec,nkin)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(OUT)                                       :: nexchange
CHARACTER (LEN=mls), INTENT(IN)                                 :: data1
INTEGER(I4B), INTENT(OUT)                                       :: nexch_sec
INTEGER(I4B), INTENT(IN)                                        :: nkin

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: dummy1
CHARACTER (LEN=mls)                                             :: namtmp

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: ltemp
INTEGER(I4B)                                                    :: lmin
INTEGER(I4B)                                                    :: k
INTEGER(I4B)                                                    :: lexc
INTEGER(I4B)                                                    :: n
INTEGER(I4B)                                                    :: ick
INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: ix
INTEGER(I4B)                                                    :: ii

REAL(DP)                                                        :: exchange_tmp
REAL(DP)                                                        :: bfit_tmp


CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                  :: nam
REAL(DP), DIMENSION(:), ALLOCATABLE                             :: sto

ALLOCATE(nam(50))
ALLOCATE(sto(50))

REWIND nout

nexchange = 0
nexch_sec = 0
iexc = 1

100 READ(nout,'(a)',END=111) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'exchange') THEN
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
!            call convan(ssch,lzs,res)
      CALL stringtype(ssch,lzs,res)
      ltemp = lzs
      IF (res == 'a') THEN
        namtmp = ssch
        nexchange = nexchange + 1
        namexc(nexchange) = namtmp
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Looking for ASCII string for name of exchanger'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF (ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch /= 'on') THEN
          WRITE(*,*) ' Name of exchanger should be followed by "on" or nothing at all'
          WRITE(*,*) ' Exchange species: ',namexc(nexchange)
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
          WRITE(*,*) ' In ion exchange block'
          WRITE(*,*) ' Exchange species: ',namexc(nexchange)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        namtmp = ssch
        DO k = 1,nkin
          IF (namtmp == umin(k)) THEN
            kexch(nexchange) = k
            GO TO 50
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,*) ' Mineral substrate listed in ion exchange block'
        WRITE(*,*) '     not found in minerals list'
        WRITE(*,*) ' Exchange species: ',namexc(nexchange)
        WRITE(*,*) ' Looking for mineral: ',namtmp(1:lmin)
        WRITE(*,*)
        READ(*,*)
        STOP
        
        50 iexchange(nexchange) = 1     ! Exchange on a specific mineral
        GO TO 100
        
      ELSE
        iexchange(nexchange) = 0     ! Exchange on bulk sediment
      END IF
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Exchange keyword specified, but no exchanger name given'
      WRITE(*,*) ' In ion exchange block'
      WRITE(*,*) ' Exchanger name must be specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  END IF
END IF
GO TO 100

!  Now, look for the activity convention to be used

111 REWIND nout

200 READ(nout,'(a)',END=222) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'convention') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      lexc = lzs
      namtmp = ssch
      IF (namtmp == 'gaines-thomas') THEN
        iexc = 1
      ELSE IF (namtmp == 'vanselow') THEN
        iexc = 2
      ELSE IF (namtmp == 'gapon') THEN
        iexc = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont understand exchanger activity convention'
        WRITE(*,*) namtmp(1:lexc)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No convention for exchanger activity given'
      WRITE(*,*) ' Assuming Gaines-Thomas convention'
      iexc = 1
      WRITE(*,*)
    END IF
  END IF
END IF
GO TO 200


222 OPEN(UNIT=18,FILE=data1,STATUS='old')

!  Find the beginning of the exchange section

300 READ(18,'(a)') dummy1
IF (dummy1 == 'Begin exchange') THEN
  GO TO 400
ELSE
  GO TO 300
END IF

400 READ(18,'(a)',END=444) dummy1
IF (dummy1 == 'End of exchange') THEN
  CLOSE(UNIT=18,STATUS='keep')
  DEALLOCATE(nam)
  DEALLOCATE(sto)
  RETURN
ELSE
  BACKSPACE 18
  READ(18,*,ERR=6003) nam(1),n,(sto(i+1),nam(i+1), i = 1, n), exchange_tmp,bfit_tmp
  
  DO ick = 1,n
    
    DO i = 1,ncomp
      IF (nam(ick+1) == ulab(i)) THEN
        GO TO 555
      END IF
    END DO
    DO ix = 1,nexchange
      IF (nam(ick+1) == namexc(ix)) THEN
        GO TO 555
      END IF
    END DO
!  Species not found in primary species or exchanger lists, skip reaction
    GO TO 400
    555     CONTINUE   ! Goto to next species in reaction
    
  END DO
  
  nexch_sec = nexch_sec + 1
  nam_exchsec(nexch_sec) = nam(1)
  DO i = 1,ncomp
    DO ick = 1,n
      IF (nam(ick+1) == ulab(i)) THEN
        muexc(nexch_sec,i) = sto(ick+1)
      END IF
    END DO
  END DO
  DO ix = 1,nexchange
    ii = ix + ncomp
    DO ick=1,n
      IF (nam(ick+1) == namexc(ix)) THEN
        ixlink(nexch_sec) = ix
        muexc(nexch_sec,ii) = sto(ick+1)
      END IF
    END DO
  END DO
  DO i = 1,ncomp
    DO ick=1,n
      IF (nam(ick+1) == ulab(i)) THEN
        nclink(nexch_sec) = i
      END IF
    END DO
  END DO
  keqexc(nexch_sec) = exchange_tmp*clg
  bfit(nexch_sec) = bfit_tmp*clg
  
END IF

GO TO 400

444 CLOSE(UNIT=18,STATUS='keep')

DEALLOCATE(nam)
DEALLOCATE(sto)

RETURN

6003 WRITE(*,*) ' Error in reading ion exchange reactions'
READ(*,*)
STOP

END SUBROUTINE read_exchange
