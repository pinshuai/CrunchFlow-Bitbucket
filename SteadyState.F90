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

    
SUBROUTINE SteadyState(ncomp,nx,ny,nz,delt,steady)
USE crunchtype
USE concentration
USE runtime

IMPLICIT NONE

!   External variables and arrays

INTEGER(I4B), INTENT(IN)                                         :: ncomp
INTEGER(I4B), INTENT(IN)                                         :: nx
INTEGER(I4B), INTENT(IN)                                         :: ny
INTEGER(I4B), INTENT(IN)                                         :: nz
REAL(DP), INTENT(IN)                                             :: delt
LOGICAL(LGT), INTENT(OUT)                                        :: steady

!  Internal variables and arrays

INTEGER(I4B)                                                     :: jx
INTEGER(I4B)                                                     :: jy
INTEGER(I4B)                                                     :: jz
INTEGER(I4B)                                                     :: i
INTEGER(I4B)                                                     :: imax
INTEGER(I4B)                                                     :: jxmax
INTEGER(I4B)                                                     :: jymax
INTEGER(I4B)                                                     :: jzmax
INTEGER(I4B)                                                     :: jxwidth
INTEGER(I4B)                                                     :: jywidth
INTEGER(I4B)                                                     :: jzwidth
INTEGER(I4B)                                                     :: width
INTEGER(I4B)                                                     :: ls

REAL(DP)                                                         :: tderivative
REAL(DP)                                                         :: maxchange

maxchange = 0.0
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO i = 1,ncomp
         tderivative = (s(i,jx,jy,jz) - sn(i,jx,jy,jz))/(delt*sn(i,jx,jy,jz))
         IF (ABS(tderivative) > maxchange) THEN
           maxchange = ABS(tderivative)
           imax = i
           jxmax = jx
           jymax = jy
           jzmax = jz
         END IF
      END DO
    END DO
  END DO
END DO

CALL stringlen(ulab(imax),ls)

WRITE(*,*) 
WRITE(*,FMT=100) maxchange
WRITE(*,*) 'Component                                    :  ', ulab(imax)(1:ls)
WRITE(*,FMT=102) jxmax,jymax,jzmax
WRITE(*,*)

100 FORMAT(' Maximum rate of change (normalized mol/kg/yr): ',1pe12.4)
102 FORMAT(' Grid location                                : ',I3,':',I3,':',I3)

IF (maxchange < steadytol) THEN
  steady = .true.
ELSE
  steady = .false.
END IF

RETURN
END SUBROUTINE SteadyState
