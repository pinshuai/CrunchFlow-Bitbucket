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

SUBROUTINE read_permz(nout,nx,ny,nz,npermz)
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
INTEGER(I4B), INTENT(OUT)                                   :: npermz

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

REAL(DP)                                                    :: permz_tmp

nxyz = nx*ny*nz

permzonez(0) = 0.0
REWIND nout

npermz = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'permeability_z') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        permz_tmp = DNUM(ssch)
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
            permzonez(0) = permz_tmp
          ELSE IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            npermz = npermz + 1
            IF (npermz > mperm) THEN
              WRITE(*,*)
              WRITE(*,*)  ' Number of permeability zones dimensioned too small'
              WRITE(*,*)  ' Number of permeability zones = ',npermz
              WRITE(*,*)  ' Dimension of permeability zones = ',mperm
              WRITE(*,*)  ' Contact the code developer at CISteefel@lbl.gov'
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            
            permzonez(npermz) = permz_tmp
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jxxpermz_lo(npermz) = JNUM(ssch_a)
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
                  jxxpermz_hi(npermz) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jxxpermz_hi(npermz) = jxxpermz_lo(npermz)   !  Assume jxxpermz_hi=jxxpermz_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for permeability'
              WRITE(*,*) ' Permeability zone ',npermz
              WRITE(*,*)
              STOP
            END IF
            
            WRITE(*,*)
            WRITE(*,*) ' Z permeability zone number ',npermz
            WRITE(*,*) ' Jxxpermz_lo = ', jxxpermz_lo(npermz)
            WRITE(*,*) ' Jxxpermz_hi = ',jxxpermz_hi(npermz)
            WRITE(*,*)

            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jyypermz_lo(npermz) = JNUM(ssch_a)
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
                  jyypermz_hi(npermz) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jyypermz_hi(npermz) = jyypermz_lo(npermz)   !  Assume jxxpermz_hi=jxxpermz_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Y location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermz
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jyypermz_lo = ',jyypermz_lo(npermz)
            WRITE(*,*) ' Jyypermz_hi = ',jyypermz_hi(npermz)
            WRITE(*,*)    

            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jzzpermz_lo(npermz) = JNUM(ssch_a)
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
                  jzzpermz_hi(npermz) = DNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jzzpermz_hi(npermz) = jzzpermz_lo(npermz)   !  Assume jxxpermz_hi=jxxpermz_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Z location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermz
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jzzpermz_lo = ',jzzpermz_lo(npermz)
            WRITE(*,*) ' Jzzpermz_hi = ',jzzpermz_hi(npermz)
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
        permzonez(0) = permz_tmp
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

500 DO l = 1,npermz
  IF (jxxpermz_hi(l) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JX > NX'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermz_hi(l) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JY > NY'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermz_hi(l) > nz+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JZ > NZ+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxpermz_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JX < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermz_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JY < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermz_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an Z permeability at JZ < 0'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_permz
