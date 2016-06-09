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


SUBROUTINE readModFlowStress(nout,perlen,nstp,tsmult,ndimdummy)
USE crunchtype
USE CrunchFunctions
USE strings
USE modflowModule

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: perlen
INTEGER(I4B), DIMENSION(:), INTENT(OUT)                     :: nstp
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: tsmult
INTEGER(I4B), INTENT(IN)                                    :: ndimdummy

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: i

REWIND nout

nstress = 0

10  READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ssch == 'modflow_stress') THEN
  nstress = nstress + 1
  IF (nstress > ndimdummy) THEN
    WRITE(*,*)
    WRITE(*,*) ' Too many stress periods'
    WRITE(*,*) ' Redimension "ndimdummy" in routine "readModFlowStress" '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Read the length of this stress period

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      perlen(nstress) = DNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading length of stress period in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading length of stress period in MODFLOW block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Now, read the number of time steps in the stress period

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      nstp(nstress) = JNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading number of time steps in stress period in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading number of time steps in stress period'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

! Now, read the multiplier

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (res == 'n') THEN
      tsmult(nstress) = DNUM(ssch)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' --> Error reading time step multiplier in MODFLOW block'
      WRITE(*,*) '     Cannot interpret string: ',ssch(1:ls)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE                  ! Zero length trailing string
    WRITE(*,*)
    WRITE(*,*) ' --> Error reading time step multiplier in MODFLOW block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

END IF

GO TO 10

500  RETURN

END SUBROUTINE readModFlowStress
