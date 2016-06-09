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


SUBROUTINE read_speciesdiffusion(nout,ncomp,nspec,ndiff)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(OUT)                                   :: ndiff

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: spectemp
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: iksave
INTEGER(I4B)                                                :: lspec


REWIND nout

ndiff = 0

10  READ(nout,'(a)',END=1000) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch /= 'd_25') THEN
    GO TO 10
  END IF
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  spectemp = ssch
  CALL stringlen(spectemp,lspec)
  IF (ls == 0) THEN
    GO TO 10
  END IF
  DO ik = 1,ncomp+nspec
    IF (spectemp == ulab(ik)) THEN
      iksave = ik
      GO TO 50
    END IF
  END DO
 
  WRITE(*,*)
  WRITE(*,*) ' Looking for diffusion coefficients for individual species'
  WRITE(*,*) ' Species not found: ',spectemp(1:lspec)
  WRITE(*,*) ' In transport block'
  WRITE(*,*)
  READ(*,*)
  STOP
  
  50    id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      ndiff = ndiff + 1
      d_sp(iksave) = DNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' Cant interpret string following "D_25 and the species name" '
      WRITE(*,*) ' A numerical value should follow'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No value for diffusion coefficient following species name'
    WRITE(*,*) ' Species = ',spectemp(1:lspec)
    WRITE(*,*) ' Assuming default value'
    WRITE(*,*)
  END IF
  
ELSE
  GO TO 10
END IF

GO TO 10

1000 RETURN
END SUBROUTINE read_speciesdiffusion
