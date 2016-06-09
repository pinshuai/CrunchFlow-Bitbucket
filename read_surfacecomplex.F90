!! CrunchTope 
!! Copyright (c) 2016, Carl Steefel
!! Copyright (c) 2016, The Regents of the University of California, 
!! through Lawrence Berkeley National Laboratory (subject to 
!! receipt of any required approvals from the U.S. Dept. of Energy).  
!! All rights reserved.

!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met: 

!! (1) Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.

!! (2) Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.

!! (3) Neither the name of the University of California, Lawrence
!! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of    
!! its contributors may be used to endorse or promote products derived
!! from this software without specific prior written permission.

!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE 

    
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
complexfound = .FALSE.

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  IF (ssch == namsurf(ks)) THEN
    complexfound = .TRUE.     ! Surface complex found
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
