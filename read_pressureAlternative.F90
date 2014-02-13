!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:06:39
 
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

SUBROUTINE read_pressureAlternative(nout,nx,ny,nz,npressure)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: npressure

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: l
INTEGER(I4B)                                                :: ArraySize

REAL(DP)                                                    :: PressTmp

ArraySize = 1000

ALLOCATE(PressureZone(0:ArraySize))
ALLOCATE(PressureFix(ArraySize))
ALLOCATE(jxxPressure_lo(ArraySize))
ALLOCATE(jxxPressure_hi(ArraySize))
ALLOCATE(jyyPressure_lo(ArraySize))
ALLOCATE(jyyPressure_hi(ArraySize))
ALLOCATE(jzzPressure_lo(ArraySize))
ALLOCATE(jzzPressure_hi(ArraySize))

PressureFix = 1
PressureZone = 0.000
jxxPressure_lo = 1
jxxPressure_hi = 1
jyyPressure_lo = 1
jyyPressure_hi = 1
jzzPressure_lo = 1
jzzPressure_hi = 1

nxyz = nx*ny*nz

npressure = 0
REWIND nout

10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'pressure') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        PressTmp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "pressure"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        STOP
      END IF
      
! Now look for ASCII string indicating location of pressure
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'a') THEN
          IF (ssch == 'default' .OR. ssch == 'all') THEN
            PressureZone(0) = PressTmp
          ELSE IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            npressure = npressure + 1
            IF (npressure > 1000) THEN
              WRITE(*,*)
              WRITE(*,*)  ' Number of pressure zones dimensioned too small'
              WRITE(*,*)  ' Number of pressure zones = ', npressure
              WRITE(*,*)  ' Dimension of pressure zones = 1000'
              WRITE(*,*)
              STOP
            END IF
            
            PressureZone(npressure) = PressTmp
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jxxPressure_lo(npressure) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string'
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jxxPressure_hi(npressure) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "pressure"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jxxPressure_hi(npressure) = jxxPressure_lo(npressure)   !  Assume jxxpermx_hi=jxxpermx_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for pressure'
              WRITE(*,*) ' Pressure zone ', npressure
              WRITE(*,*)
              STOP
            END IF
            
!!            IF (ny > 1) THEN
              id = ids + ls
              CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
              IF(ls /= 0) THEN
                lzs=ls_a
                CALL convan(ssch_a,lzs,res)
                IF (res == 'n') THEN
                  jyyPressure_lo(npressure) = JNUM(ssch_a)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' No Y location for pressure '
                  WRITE(*,*)
                  STOP
                END IF
                IF (ls_b /= 0) THEN
                  lzs=ls_b
                  CALL convan(ssch_b,lzs,res)
                  IF (res == 'n') THEN
                    jyyPressure_hi(npressure) = JNUM(ssch_b)
                  ELSE                !  An ascii string--so bag it.
                    WRITE(*,*)
                    WRITE(*,*) ' A grid location should follow zone specification'
                    WRITE(*,*) ' Dont know what to do with this string after "pressure"'
                    WRITE(*,*)
                    STOP
                  END IF
                ELSE
                  jyyPressure_hi(npressure) = jyyPressure_lo(npressure)   !  Assume jxxpermx_hi=jxxpermx_lo
                END IF
              ELSE                  ! Zero length trailing string
                WRITE(*,*)
                WRITE(*,*) ' No Y location given for pressure zone'
                WRITE(*,*) ' Pressure zone number ', npressure
                WRITE(*,*)
                STOP
              END IF
              
              
!!            ELSE
!!              jyyPressure_lo(npressure) = 1
!!              jyyPressure_hi(npressure) = 1
!!            END IF

!!            IF (nz > 1) THEN
              id = ids + ls
              CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
              IF(ls /= 0) THEN
                lzs=ls_a
                CALL convan(ssch_a,lzs,res)
                IF (res == 'n') THEN
                  jzzPressure_lo(npressure) = JNUM(ssch_a)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' No Z location for pressure '
                  WRITE(*,*)
                  STOP
                END IF
                IF (ls_b /= 0) THEN
                  lzs=ls_b
                  CALL convan(ssch_b,lzs,res)
                  IF (res == 'n') THEN
                    jzzPressure_hi(npressure) = JNUM(ssch_b)
                  ELSE                !  An ascii string--so bag it.
                    WRITE(*,*)
                    WRITE(*,*) ' A grid location should follow zone specification'
                    WRITE(*,*) ' Dont know what to do with this string after "pressure"'
                    WRITE(*,*)
                    STOP
                  END IF
                ELSE
                  jzzPressure_hi(npressure) = jzzPressure_lo(npressure)   !  Assume jxxpermx_hi=jxxpermx_lo
                END IF
              ELSE                  ! Zero length trailing string
                WRITE(*,*)
                WRITE(*,*) ' No Z location given for pressure zone'
                WRITE(*,*) ' Pressure zone number ', npressure
                WRITE(*,*)
                STOP
              END IF          
              
!!            ELSE
!!              jzzPressure_lo(npressure) = 1
!!              jzzPressure_hi(npressure) = 1
!!            END IF  
            
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following pressure value'
            WRITE(*,*) ssch(1:ls)
            WRITE(*,*)
            STOP
          END IF

          id = ids + ls
          CALL sschaine(zone,id,iff,ssch,ids,ls)
          IF (ls /= 0) THEN
            lzs=ls
            CALL convan(ssch,lzs,res)
            IF (ssch == 'fix') THEN
              PressureFix(npressure) = 0
            ELSE
              PressureFix(npressure) = 1
            END IF
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following pressure value'
          WRITE(*,*) ' Looking for ASCII string'
          WRITE(*,*)
          STOP
        END IF
      ELSE   ! Assume this is default if nothing else given
        PressureZone(0) = PressTmp
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value given for pressure'
      WRITE(*,*) ' Pressure specification ignored'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    GO TO 10
  END IF
  
END IF

GO TO 10

500 DO l = 1,npressure
  IF (jxxPressure_hi(l) > nx+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JX > nx+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyPressure_hi(l) > ny+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JY > ny+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzPressure_hi(l) > nz+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JZ > nz+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxPressure_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JX < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyPressure_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JY < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzPressure_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an pressure at JZ < 0'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_pressureAlternative
