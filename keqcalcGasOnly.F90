!******************        GIMRT98      ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:00:11
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!********************** C.I. Steefel *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE keqcalcGasOnly(ngas,nspec,jx,jy,jz,pg)
USE crunchtype
USE params
USE concentration
USE mineral
USE temperature

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ngas
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: jx
INTEGER(I4B), INTENT(IN)                             :: jy
INTEGER(I4B), INTENT(IN)                             :: jz
REAL(DP), INTENT(IN)                                 :: pg

!  Internal variables

REAL(DP)                                             :: temp
REAL(DP)                                             :: temp2
REAL(DP)                                             :: x1
REAL(DP)                                             :: x2
REAL(DP)                                             :: x3
REAL(DP)                                             :: x4
REAL(DP)                                             :: x5
REAL(DP)                                             :: mu

INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: kg
INTEGER(I4B)                                         :: k
INTEGER(I4B)                                         :: msub
INTEGER(I4B)                                         :: ns
INTEGER(I4B)                                         :: np


temp = t(jx,jy,jz) + 273.15
temp2 = temp*temp

DO kg = 1,ngas
  ksp = kg + nspec

if_co2: if (namg(kg) == 'CO2(g)') then

  call calc_mu(pg,temp,mu)
  keqgas(kg,jx,jy,jz) = mu
!!  write(*,*)
!!  write(*,*)'updated CO2(g) equilibrium constant (log10K): ',keqgas_tmp(kg)/clg 
!!  write(*,*)
else
  
  IF (ntemp == 1 .OR. RunIsothermal) THEN
    keqgas(kg,jx,jy,jz) = -clg*eqgas(kg)
  ELSE
    x1 = as1(ksp,1)
    x2 = as1(ksp,2)
    x3 = as1(ksp,3)
    x4 = as1(ksp,4)
    x5 = as1(ksp,5)
    keqgas(kg,jx,jy,jz) = -clg*(x1*DLOG(temp) + x2 +  &
        x3*temp + x4/temp + x5/(temp2))
  END IF
  
end if if_co2

end do

RETURN
END SUBROUTINE keqcalcGasOnly
!*********************************************************************
