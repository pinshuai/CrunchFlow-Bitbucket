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
    
SUBROUTINE surf_init(ncomp,nspec,nsurf,nsurf_sec,nco)
USE crunchtype
USE concentration
USE mineral

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: nco


!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents
REAL(DP)                                                    :: CheckWrite

DO ns = 1,nsurf_sec

  IF (nptlink(ns) /= 0) THEN

    delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

!!  Aqueous species
    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
      sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
      sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))
    IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.0d0
    ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.693147180559945
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
      LogTotalEquivalents = LogTotalSites - 1.09861228866811
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
      LogTotalEquivalents = 0.0d0
    ELSE
      write(*,*)
      write(*,*) ' Do not recognize denticity in surface complexation'
      write(*,*)
      read(*,*)
      stop
    END IF

!!  Surface complexes
    DO is = 1,nsurf
      activity = spsurftmp(is) - LogTotalSites
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO

    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) -              &
       delta_z*2.0*LogPotential_tmp(nptlink(ns)) + sum + LogTotalEquivalents
    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))

  ELSE                                  !! Non-electrostatic option

!!  Aqueous species
    sum = 0.0
    DO i = 1,ncomp
      IF (ulab(i) == 'H2O') THEN
        sum = sum + musurf(ns,i)*(gamtmp(i))
      ELSE
        sum = sum + musurf(ns,i)*(sptmp(i) + gamtmp(i))
      END IF
    END DO

    LogTotalSites = DLOG(c_surf(islink(ns),nco))
    IF (musurf(ns,islink(ns)+ncomp) == 1.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.0d0
    ELSe IF (musurf(ns,islink(ns)+ncomp) == 2.0) THEN
      LogTotalEquivalents = LogTotalSites - 0.693147180559945
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 3.0) THEN
      LogTotalEquivalents = LogTotalSites - 1.09861228866811
    ELSE IF (musurf(ns,islink(ns)+ncomp) == 0.0) THEN
      LogTotalEquivalents = 0.0d0
    ELSE
      write(*,*)
      write(*,*) ' Do not recognize denticity in surface complexation'
      write(*,*)
      read(*,*)
      stop
    END IF

!!  Surface complexes
    DO is = 1,nsurf
      activity = spsurftmp(is) - LogTotalSites
      sum = sum + musurf(ns,is+ncomp)*activity
    END DO

    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum + LogTotalEquivalents
!!    spsurftmp(ns+nsurf) = keqsurf_tmp(ns) + sum + LogTotalSites
    spsurftmp10(ns+nsurf) = DEXP(spsurftmp(ns+nsurf))
  END IF

END DO

RETURN
END SUBROUTINE surf_init
!  **************************************************************
