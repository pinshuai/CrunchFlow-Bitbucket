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


SUBROUTINE porperm(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE flow
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nx
INTEGER(I4B), INTENT(IN)                                        :: ny
INTEGER(I4B), INTENT(IN)                                        :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                    :: jx
INTEGER(I4B)                                                    :: jy
INTEGER(I4B)                                                    :: jz

permxOld = permx
permyOld = permy
permzOld = permz

do jz = 1,nz
DO jy = 1,ny
  DO jx = 1,nx
!!    permx(jx,jy,jz) = permxOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
!!                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 ) 
    permx(jx,jy,jz) = permxOld(jx,jy,jz)*( por(jx,jy,jz)/porOld(jx,jy,jz) )**3  * ( (1.0d0-porOld(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0
    permy(jx,jy,jz) = permyOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 )  
    permz(jx,jy,jz) = permzOld(jx,jy,jz) * ( por(jx,jy,jz)          )**3.0/ ( ( 1.0d0-por(jx,jy,jz) )**2.0 ) * &
                                           ( 1.0d0-porOld(jx,jy,jz) )**2.0/ ( (porOld(jx,jy,jz)     )**3.0 )  
    if (jx==12) then
        continue
    end if
!!    permx(jx,jy,jz) = perminx(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0
!!    permy(jx,jy,jz) = perminy(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0 
!!    permz(jx,jy,jz) = perminz(jx,jy,jz)*( por(jx,jy,jz)/porin(jx,jy,jz) )**3  * ( (1.0d0-porin(jx,jy,jz))/(1.0d0-por(jx,jy,jz)) )**2.0d0 
  END DO
END DO
end do


!!DO jy = 1,ny
!!  DO jx = 1,nx
!!    IF (por(jx,jy,jz) > 0.20) THEN
!!      permx(jx,jy,jz) = perminx(jx,jy,jz)*10000000.0
!!      permy(jx,jy,jz) = perminy(jx,jy,jz)*10000000.0
!!      permz(jx,jy,jz) = perminz(jx,jy,jz)*10000000.0
!!    ELSE
!!      permx(jx,jy,jz) = perminx(jx,jy,jz)
!!      permy(jx,jy,jz) = perminy(jx,jy,jz)
!!      permz(jx,jy,jz) = perminz(jx,jy,jz)
!!    END IF!!
!!  END DO
!!END DO

RETURN
END SUBROUTINE porperm

