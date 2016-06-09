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
    
SUBROUTINE SurfLocal(ncomp,nsurf,nsurf_sec,jx,jy,jz,AqueousToBulk)
USE crunchtype
USE concentration
USE mineral
USE medium
USE transport
USE temperature
USE RunTime

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                    :: jx
INTEGER(I4B), INTENT(IN)                                    :: jy
INTEGER(I4B), INTENT(IN)                                    :: jz
REAL(DP), INTENT(IN)                                        :: AqueousToBulk

!  Internal variables and arrays

INTEGER(I4B)                                                :: ns
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: is

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: delta_z
REAL(DP)                                                    :: LogAqueousToBulk
REAL(DP)                                                    :: activity
REAL(DP)                                                    :: LogTotalSites
REAL(DP)                                                    :: LogTotalEquivalents


!!      CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)
!!      LogAqueousToBulk = DLOG(AqueousToBulk)

      DO ns = 1,nsurf_sec

        IF (nptlink(ns) /= 0) THEN                            !  Electrostatic correction

          delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF

          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 
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

          DO is = 1,nsurf

!!            activity = spsurf(is,jx,jy,jz) - LogTotalSites - LogAqueousToBulk
            activity = spsurf(is,jx,jy,jz) - LogTotalSites
            sum = sum + musurf(ns,is+ncomp)*activity

          END DO

!! NOTE:  Below is the LOg concentration of sites in units of mol/kgw
          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum -                     &
             delta_z*2.0*LogPotential(nptlink(ns),jx,jy,jz) + LogTotalEquivalents

!!          spsurf10(ns+nsurf,jx,jy,jz) = AqueousToBulk*DEXP( spsurf(ns+nsurf,jx,jy,jz) )
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        ELSE                                                  !  Non-electrostatic 

          IF (ikh2o /= 0) THEN

            sum = 0.0
            DO i = 1,ncomp
              IF (ulab(i) == 'H2O') THEN
                sum = sum + musurf(ns,i)*(gam(i,jx,jy,jz))
              ELSE
                sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
              END IF
            END DO

          ELSE

            sum = 0.0
            DO i = 1,ncomp
              sum = sum + musurf(ns,i)*(sp(i,jx,jy,jz) + gam(i,jx,jy,jz))
            END DO

          END IF

          LogTotalSites = LogTotalSurface(islink(ns),jx,jy,jz) 
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

          DO is = 1,nsurf

            activity = spsurf(is,jx,jy,jz) - LogTotalSites
            sum = sum + musurf(ns,is+ncomp)*activity

          END DO

          spsurf(ns+nsurf,jx,jy,jz) = keqsurf(ns,jx,jy,jz) + sum + LogTotalEquivalents

!!          spsurf10(ns+nsurf,jx,jy,jz) = AqueousToBulk*DEXP( spsurf(ns+nsurf,jx,jy,jz) )
          spsurf10(ns+nsurf,jx,jy,jz) = DEXP( spsurf(ns+nsurf,jx,jy,jz) )

        END IF
      END DO
RETURN
END SUBROUTINE SurfLocal
!  **************************************************************
