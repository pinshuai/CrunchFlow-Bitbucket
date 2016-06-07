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

SUBROUTINE totconc_plus(ncomp,nspec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE transport

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz

!  Internal variables and arrays

INTEGER(I4B)                                                :: ksp
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: kk

REAL(DP)                                                    :: sum1
REAL(DP)                                                    :: sum2
REAL(DP)                                                    :: sum3
REAL(DP)                                                    :: sum4
REAL(DP)                                                    :: tmp

sum4 = 0.0
DO i = 1,ncomp
  sum1 = 0.0
  sum2 = 0.0
  sum3 = 0.0
  DO ksp = 1,nspec
    kk = ksp + ncomp
    tmp = muaq(ksp,i)*sp10(kk,jx,jy,jz)
    sum1 = sum1 + tmp
    sum2 = sum2 + tmp*d_correct(kk)
    sum3 = sum3 + tmp*d_correct(kk)*chg(kk)
  END DO
  tmp = sp10(i,jx,jy,jz)
  s(i,jx,jy,jz)     = sum1 + tmp
  s_dsp(i,jx,jy,jz) = sum2 + tmp*d_correct(i)
  s_chg(i,jx,jy,jz) = sum3 + tmp*d_correct(i)*chg(i)
  sum4 = sum4 + s_chg(i,jx,jy,jz)*chg(i)
END DO

sumwtchg(jx,jy,jz) = sum4

RETURN
END SUBROUTINE totconc_plus
!  **************************************************************
