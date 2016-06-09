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


SUBROUTINE porcalc(nx,ny,nz,nkin,jpor)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz
INTEGER(I4B), INTENT(IN)                                        :: nkin
INTEGER(I4B), INTENT(IN)                                        :: jpor

!  Internal variables and arrays

REAL(DP)                                                        :: sum
REAL(DP)                                                        :: vinit
REAL(DP)                                                        :: porfactor

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz
INTEGER(I4B)                                                    :: k

!  NOTE:  Initial volume fractions should be advected as well
!         Assumed immobile for the moment

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      sum = 0.0
      DO k = 1,nkin
        vinit = volin(k,jinit(jx,jy,jz))
        IF (vinit == 0.0) THEN
          area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))*(volfx(k,jx,jy,jz)/0.01)**0.66666
          if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz)
        ELSE
          if (mintype(k) == 0) sum = sum + volfx(k,jx,jy,jz)
          area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))* (volfx(k,jx,jy,jz)/vinit)**0.6666
        END IF
      END DO
      IF (jpor /= 1) THEN
        por(jx,jy,jz) = porin(jx,jy,jz)
      ELSE
        por(jx,jy,jz) = 1.0 - sum
        IF (por(jx,jy,jz) <= 0.0) THEN
          por(jx,jy,jz) = 1.d-14
        END IF
      END IF
      porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**0.6666666666666
      DO k = 1,nkin
        area(k,jx,jy,jz) = area(k,jx,jy,jz)*porfactor
      END DO
    END DO
  END DO
END DO

RETURN
END SUBROUTINE porcalc
!  *******************************************************
