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

SUBROUTINE read_permy(nout,nx,ny,nz,npermy)
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
INTEGER(I4B), INTENT(OUT)                                   :: npermy

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

REAL(DP)                                                    :: permy_tmp

nxyz = nx*ny*nz

permzoney(0) = 0.0
REWIND nout

npermy = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'permeability_y') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        permy_tmp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "permeability"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        STOP
      END IF
      
! Now look for ASCII string indicating location of permeability
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'a') THEN
          IF (ssch == 'default' .OR. ssch == 'all') THEN
            permzoney(0) = permy_tmp
          ELSE IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            npermy = npermy + 1
            IF (npermy > mperm) THEN
              WRITE(*,*)
              WRITE(*,*)  ' Number of permeability zones dimensioned too small'
              WRITE(*,*)  ' Number of permeability zones = ',npermy
              WRITE(*,*)  ' Dimension of permeability zones = ',mperm
              WRITE(*,*)  ' Contact the code developer at CISteefel@lbl.gov'
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            
            permzoney(npermy) = permy_tmp
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
               jxxpermy_lo(npermy) = JNUM(ssch_a)
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
                  jxxpermy_hi(npermy) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jxxpermy_hi(npermy) = jxxpermy_lo(npermy)   !  Assume jxxpermy_hi=jxxpermy_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for permeability'
              WRITE(*,*) ' Permeability zone ',npermy
              WRITE(*,*)
              STOP
            END IF
            
            WRITE(*,*)
            WRITE(*,*) ' Y permeability zone number ',npermy
            WRITE(*,*) ' Jxxpermy_lo = ', jxxpermy_lo(npermy)
            WRITE(*,*) ' Jxxpermy_hi = ',jxxpermy_hi(npermy)
            WRITE(*,*)
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jyypermy_lo(npermy) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Y location for permeability '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jyypermy_hi(npermy) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jyypermy_hi(npermy) = jyypermy_lo(npermy)   !  Assume jxxpermy_hi=jxxpermy_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Y location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermy
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jyypermy_lo = ',jyypermy_lo(npermy)
            WRITE(*,*) ' Jyypermy_hi = ',jyypermy_hi(npermy)
            WRITE(*,*)    
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jzzpermy_lo(npermy) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Z location for permeability '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jzzpermy_hi(npermy) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jzzpermy_hi(npermy) = jzzpermy_lo(npermy)   !  Assume jxxpermy_hi=jxxpermy_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Z location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermy
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jzzpermy_lo = ',jzzpermy_lo(npermy)
            WRITE(*,*) ' Jzzpermy_hi = ',jzzpermy_hi(npermy)
            WRITE(*,*)    

          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following permeability value'
            WRITE(*,*) ssch(1:ls)
            WRITE(*,*)
            STOP
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following permeability value'
          WRITE(*,*) ' Looking for ASCII string'
          WRITE(*,*)
          STOP
        END IF
      ELSE   ! Assume this is default if nothing else given
        permzoney(0) = permy_tmp
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value given for permeability'
      WRITE(*,*) ' Permeability specification ignored'
      WRITE(*,*)
    END IF
  ELSE
    GO TO 10
  END IF
  
END IF

GO TO 10

500 DO l = 1,npermy
  IF (jxxpermy_hi(l) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JX > NX'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermy_hi(l) > ny+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JY > NY+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermy_hi(l) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JZ > NZ'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxpermy_lo(l) < 01) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JX < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermy_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JY < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermy_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Y permeability at JZ < 1'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_permy
