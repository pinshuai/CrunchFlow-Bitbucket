!************** (C) COPYRIGHT 1993 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:43:21
 
!                      All Rights Reserved

!  GIMRT IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR
!  IMPLIED. THE USER ASSUMES ALL RISKS OF USING 1DREACT. THERE  IS
!  NO CLAIM OF THE MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

SUBROUTINE GraphicsKaleidagraph(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
    ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,delt,jpor)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE strings
USE NanoCrystal

IMPLICIT NONE
!  *********************  INTERFACE BLOCKS  *****************************
INTERFACE
  SUBROUTINE GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
    USE crunchtype
    INTEGER(I4B), INTENT(IN)                                   :: ncomp
    INTEGER(I4B), INTENT(IN)                                   :: ngas
    REAL(DP), DIMENSION(:)                                     :: gastmp10
    INTEGER(I4B), INTENT(IN)                                   :: jx
    INTEGER(I4B), INTENT(IN)                                   :: jy
    INTEGER(I4B), INTENT(IN)                                   :: jz
  END SUBROUTINE GasPartialPressure
END INTERFACE
!  **********************************************************************

!  External variables and arrays

REAL(DP), INTENT(IN)                               :: realtime
REAL(DP), INTENT(IN)                               :: delt

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: ndecay
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: ikin
INTEGER(I4B), INTENT(IN)                           :: nkin
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz
INTEGER(I4B), INTENT(IN)                           :: nn
INTEGER(I4B), INTENT(IN)                           :: nint
INTEGER(I4B), INTENT(IN)                           :: ikmast
INTEGER(I4B), INTENT(IN)                           :: ikph
INTEGER(I4B), INTENT(IN)                           :: jpor

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=mls)                                 :: fn
CHARACTER (LEN=mls)                                  :: suf
CHARACTER (LEN=mls)                                  :: suf1
CHARACTER (LEN=mls)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
 
INTEGER(I4B), DIMENSION(nrct)                      :: len_min
INTEGER(I4B)                                       :: j
INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz
INTEGER(I4B)                                       :: ilength
INTEGER(I4B)                                       :: ik
INTEGER(I4B)                                       :: k
INTEGER(I4B)                                       :: ks
INTEGER(I4B)                                       :: ns
INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ir
INTEGER(I4B)                                       :: lsjx
INTEGER(I4B)                                       :: nlen
INTEGER(I4B)                                       :: l

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(ncomp)                         :: WeightPercent
REAL(DP), DIMENSION(nrct)                          :: MineralPercent
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc

REAL(DP)                                                   :: sumiap
REAL(DP)                                                   :: pHprint
REAL(DP)                                                   :: peprint
REAL(DP)                                                   :: Ehprint
REAL(DP)                                                   :: spprint
REAL(DP)                                                   :: totcharge
REAL(DP)                                                   :: siprnt
REAL(DP)                                                   :: actprint
REAL(DP)                                                   :: actprint10
REAL(DP)                                                   :: spbase
REAL(DP)                                                   :: rone
REAL(DP)                                                   :: PrintTime
REAL(DP)                                                   :: alk
REAL(DP)                                                   :: tflux_top
REAL(DP)                                                   :: tflux_bot
REAL(DP)                                                   :: top_norm
REAL(DP)                                                   :: bot_norm
REAL(DP)                                                   :: aflux_net
REAL(DP)                                                   :: ad_net_bot
REAL(DP)                                                   :: SolidRatio
REAL(DP)                                                   :: FluidRatio
REAL(DP)                                                   :: Solid234U
REAL(DP)                                                   :: Solid238U
REAL(DP)                                                   :: SolidCa
REAL(DP)                                                   :: SolidRatioMarine1
REAL(DP)                                                   :: SolidRatioMarine2
REAL(DP)                                                   :: SolidSolutionRatioTemp

REAL(DP)                                                   :: AreaWrite3
REAL(DP)                                                   :: AreaWrite4

REAL(DP)                                                   :: Del34S_sulfate
REAL(DP)                                                   :: Del34S_sulfide
REAL(DP)                                                   :: DelCa44

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: kk

REAL(DP), DIMENSION(ngas)                                  :: gastmp10
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: tk

REAL(DP)                                                   :: MoleFraction40Mineral, MoleFraction44Mineral, MoleFraction40, MoleFraction44

REAL(DP)                                                        :: pi

pi = DACOS(-1.0d0)

PrintTime = realtime*OutputTimeScale

rone = 1.0d0

suf='.out'
suf1 ='.out'
tab = CHAR(9)

DO k = 1,nrct
  uminprnt(k) = umin(k)
END DO
DO ik = 1,ncomp+nspec
  ulabprnt(ik) = ulab(ik)
END DO
DO ks = 1,nsurf
  prtsurf(ks) = namsurf(ks)
END DO
DO ns = 1,nsurf_sec
  prtsurf(ns+nsurf) = namsurf_sec(ns)
END DO

!  Write out master variable

IF (ikph /= 0) THEN
  fn='pH'
  ilength = 2
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,2282)
  jy = 1
  jz = 1
  DO jx = 1,nx
!fp! if_onproc({#expr# sp(ikph,jx,jy,jz) #});
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,183) x(jx)*OutputDistanceScale,phprt
!fp! end_onproc();
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

fn='conc'
ilength = 4
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,101)
101 FORMAT('# Units: Log10 mol/kgw')

IF (ikph /= 0) THEN
  WRITE(8,2288) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,184) x(jx)*OutputDistanceScale, phprt,(sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
  END DO
ELSE
  WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    WRITE(8,184) x(jx)*OutputDistanceScale, (sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
  END DO
END IF
CLOSE(UNIT=8,STATUS='keep')

fn='totcon'
ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,102)
102 FORMAT('# Units: Mol/kgw')
WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# s(i,jx,jy,jz) #});
  WRITE(8,184) x(jx)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

!!!fn='isotopes'
!!!ilength = 8
!!!CALL newfile(fn,suf1,fnv,nint,ilength)
!!!OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
!!!WRITE(8,2283) PrintTime
!!!WRITE(8,102)
!!!WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp)
!!!jy = 1
!!!jz = 1
!!!DO jx = 1,nx
!!!  del34S_sulfate = ( ( s(11,jx,jy,jz)/s(10,jx,jy,jz) )/0.0441626 - 1.0d0)*1000.0d0
!!!  del34S_sulfide = ( ( s(13,jx,jy,jz)/s(12,jx,jy,jz) )/0.0441626 - 1.0d0)*1000.0d0
!!!  WRITE(8,184) x(jx)*OutputDistanceScale,del34S_sulfate,del34S_sulfide
!!!END DO

!!!DO jx = 1,nx
!!!  delCa44 = ( ( s(6,jx,jy,jz)/s(7,jx,jy,jz) )/47.153 - 1.0d0)*1000.0d0
!!  del34S_sulfide = ( ( s(13,jx,jy,jz)/s(12,jx,jy,jz) )/0.0441626 - 1.0d0)*1000.0d0
!!!  WRITE(8,184) x(jx)*OutputDistanceScale,delCa44
!!!END DO

CLOSE(UNIT=8,STATUS='keep')

fn='CaIsotopes'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,102)
WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp)
jy = 1
jz = 1

!!DO jx = 1,nx
!!  delCa44 = ( ( s(6,jx,jy,jz)/s(7,jx,jy,jz) )/47.153 - 1.0d0)*1000.0d0
!!  WRITE(8,184) x(jx)*OutputDistanceScale,delCa44
!!END DO

CLOSE(UNIT=8,STATUS='keep')

fn='gas'
ilength = 3
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,103)
103 FORMAT('# Units: Bars')
WRITE(8,2285) (namg(kk),kk=1,ngas)
jy = 1
jz = 1
DO jx = 1,nx
  tk = 273.15d0 + t(jx,jy,jz)
  denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
  CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
  WRITE(8,184) x(jx)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
END DO
CLOSE(UNIT=8,STATUS='keep')

IF (nexchange > 0) THEN
  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
!!  104 FORMAT('# Units: Mol/m^3 porous medium')
  104 FORMAT('# Units: Mol/g solid')
  WRITE(8,2285) (nam_exchsec(nex),nex=1,nexch_sec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
    WRITE(8,184) x(jx)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz)/SolidSolutionRatioTemp,nex = 1,nexch_sec)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
  WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
    totex_bas = 0.0
    DO i = 1,ncomp  
      DO nex = 1,nexch_sec
        totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
      END DO
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF            !!  End of exchange block

IF (nsurf > 0) THEN
fn='surface'
ilength = 7
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,104)
WRITE(8,2285) (prtsurf(ks),ks=1,nsurf+nsurf_sec)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
  WRITE(8,184) x(jx)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz),ns = 1,nsurf+nsurf_sec)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='totsurface'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,104)
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
  totex_bas = 0.0
  DO i = 1,ncomp  
    DO ns = 1,nsurf_sec
      totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
    END DO
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nrct > 0) THEN
fn='TotMineral'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,105)
105 FORMAT('# Units: Mole component in mineral/m^3 porous medium')
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
  totex_bas = 0.0
  DO i = 1,ncomp  
    DO k = 1,nrct
      IF (volmol(k) /= 0.0) THEN
        IF (nradmax > 0) THEN
          totex_bas(i) = totex_bas(i) + mumin_decay(1,k,i,jx,1,1)*volfx(k,jx,jy,jz)/volmol(k)
        ELSE 
          totex_bas(i) = totex_bas(i) + mumin(1,k,i)*volfx(k,jx,jy,jz)/volmol(k)
        END IF
      ENDIF
    END DO
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)

END DO
CLOSE(UNIT=8,STATUS='keep')

fn='WeightPercent'
ilength = 13
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,106)
106 FORMAT('# Units: Weight % Component')
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
  totex_bas = 0.0
  sum = 0.0
  DO k = 1,nrct
    sum = sum + wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)
  END DO
  DO i = 1,ncomp  
    DO k = 1,nrct
      IF (volmol(k) /= 0.0) THEN
        IF (nradmax > 0) THEN
          totex_bas(i) = totex_bas(i) + mumin_decay(1,k,i,jx,1,1)*volfx(k,jx,jy,jz)/volmol(k)
        ELSE 
          totex_bas(i) = totex_bas(i) + mumin(1,k,i)*volfx(k,jx,jy,jz)/volmol(k)
        END IF
      ENDIF
    END DO
  END DO
  DO i = 1,ncomp
    WeightPercent(i) = 100.0*totex_bas(i)*wtcomp(i)/sum
  END DO
   
  WRITE(8,184) x(jx)*OutputDistanceScale,(WeightPercent(i),i = 1,ncomp)
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='MineralPercent'
ilength = 14
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,107)
107 FORMAT('# Units: Weight % Mineral')
WRITE(8,2285)  (uminprnt(k),k=1,nrct)
jy = 1
jz = 1
DO jx = 1,nx
  sum = 0.0
  DO k = 1,nrct
    sum = sum + wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)
  END DO 
  DO k = 1,nrct
    IF (volmol(k) /= 0.0) THEN
      MineralPercent(k) = 100.0*wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)/sum
    ENDIF
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(MineralPercent(k),k = 1,nrct)
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the reaction rates in units of mol/m**3(bulk vol.)/sec

fn='rate'
ilength = 4
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,108)
108 FORMAT('# Units: Mol/m**3 Porous Medium/s')
IF (JennyDruhan) THEN
  WRITE(8,2299) (uminprnt(k),k=1,nrct), ' 40CaBackward ', ' 40CaForward '
ELSE
  WRITE(8,2285)  (uminprnt(k),k=1,nrct)
END IF
jy = 1
jz = 1
DO jx = 1,nx
  sum = 0.0
  DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!          dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
!!        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!!        dptprt(k) = dppt(k,jx,jy,jz)/(por(jx,jy,jz))    ! mol/L fluid/yr
        dptprt(k) = dppt(k,jx,jy,jz)/(secyr)    ! mol/m**3 porous medium/sec
!*************************************
    sum = sum + dptprt(k)
  END DO
  porcalc = sum
  IF (JennyDruhan) THEN
    WRITE(8,184) x(jx)*OutputDistanceScale,(dptprt(k),k=1,nrct),rminSaveForDepaolo(1,jx,jy,jz)/secyr,rminSaveForDepaolo(2,jx,jy,jz)/secyr
  ELSE
    WRITE(8,184) x(jx)*OutputDistanceScale,(dptprt(k),k=1,nrct)
  END IF
END DO
CLOSE(UNIT=8,STATUS='keep')

ENDIF

!   Write out the reaction rates in units of mol/kgw/sec

IF (ikin > 0) THEN
fn='AqRate'
ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,109)
109 FORMAT('# Units: Mol/kgw/yr')
WRITE(8,2285)  (namkin(ir),ir=1,ikin)
jy = 1
jz = 1
DO jx = 1,nx
  sum = 0.0
  WRITE(8,184) x(jx)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
END DO
CLOSE(UNIT=8,STATUS='keep')
END IF

IF (ikin > 0) THEN
fn='DelGBiomass'
ilength = 11
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,109)
WRITE(8,2285)  (namkin(ir),ir=1,ikin)
jy = 1
jz = 1
DO jx = 1,nx
  Tk = 273.15d0 + t(jx,jy,jz)
  sum = 0.0
  WRITE(8,184) x(jx)*OutputDistanceScale,(rgas*Tk*satlog(ir,jx,jy,jz),ir=1,ikin)
END DO
CLOSE(UNIT=8,STATUS='keep')
END IF

IF (ikin > 0) THEN
fn='fTBiomass'
ilength = 9
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,109)
WRITE(8,2285)  (namkin(ir),ir=1,ikin)
jy = 1
jz = 1
DO jx = 1,nx
  Tk = 273.15d0 + t(jx,jy,jz)
  sum = 0.0
  WRITE(8,184) x(jx)*OutputDistanceScale,(1.0-DEXP(satlog(ir,jx,jy,jz)),ir=1,ikin)
END DO
CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nrct > 0) THEN

IF (JennyDruhan) then
!!  Mole fractions Ca40 and Ca44
  fn='CaIsotopeMoleFractions'
  ilength = 22
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  !!WRITE(8,*) "  Distance    Ca40FractionLT Ca44FractionLT  Ca40FractionAqueous  Ca44FractionAqueous "
  jy = 1
  jz = 1
  DO jx = 1,nx
    MoleFraction44 = sn(7,jx,jy,jz)/( sn(6,jx,jy,jz) + sn(7,jx,jy,jz) )   !! NOTE: Using old time step total concentrations (lagged)
    MoleFraction40 = 1.0d0 - MoleFraction44
    MoleFraction44Mineral = VolumeLastTimeStep(6,jx,jy,jz)/( VolumeLastTimeStep(6,jx,jy,jz) + VolumeLastTimeStep(5,jx,jy,jz) )
    MoleFraction40Mineral = VolumeLastTimeStep(5,jx,jy,jz)/( VolumeLastTimeStep(6,jx,jy,jz) + VolumeLastTimeStep(5,jx,jy,jz) )
    WRITE(8,184) x(jx)*OutputDistanceScale,MoleFraction40Mineral, MoleFraction44Mineral, MoleFraction40, MoleFraction44

  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!!  Volumes in %
  fn='volume'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,110)
  110 FORMAT('# Units: Volume % (m^3 mineral/m^3 porous medium)')
  WRITE(8,2286)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    sum = 0.0
    DO k = 1,nrct
      dvolpr(k) = 100*volfx(k,jx,jy,jz)
      sum = sum + volfx(k,jx,jy,jz)
    END DO
    porprt = (1.0-sum)*100.0
    WRITE(8,184) x(jx)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!!  Bulk areas in m2/m3 
  fn='area'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,111)
  111 FORMAT('# Units: M^2 mineral/m^3 porous medium')
  WRITE(8,2286)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    WRITE(8,184) x(jx)*OutputDistanceScale,(area(k,jx,jy,jz),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!!  Bulk areas in m2/m3 
!!CSD  fn='CSDarea'
!!CSD  ilength = 7
!!CSD  CALL newfile(fn,suf1,fnv,nint,ilength)
!!CSD  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
!!CSD  WRITE(8,2283) PrintTime
!!CSD  WRITE(8,111)
!!CSD  WRITE(8,2286)  ' Kaolinite1', 'Kaolinite2'
!!CSD  jy = 1
!!CSD  jz = 1
!!CSD  DO jx = 1,nx
!!CSD    areaWrite3 = 0.0d0
!!CSD    areaWrite4 = 0.0d0
!!CSD    DO l = 1,nCSD
!!CSD      AreaWrite3 = AreaWrite3 + NucleationScaleFactor*nCrystal(l,3,jx,jy,jz) * 2.0d0*pi*radius(l)*radius(l)
!!CSD      AreaWrite4 = AreaWrite4 + NucleationScaleFactor*nCrystal(l,4,jx,jy,jz) * 2.0d0*pi*radius(l)*radius(l)
!!CSD    END DO
!!CSD    WRITE(8,184) x(jx)*OutputDistanceScale,AreaWrite3, AreaWrite4
!!CSD  END DO
!!CSD  CLOSE(UNIT=8,STATUS='keep')

!!CSD  fn='CSD3'
!!CSD  ilength = 4
!!CSD  CALL newfile(fn,suf1,fnv,nint,ilength)
!!CSD  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
!!CSD  WRITE(8,2283) PrintTime
!!CSD  WRITE(8,111)
!!CSD  WRITE(8,*)  ' Radius    nCrystals'

!!CSD    DO l = 1,nCSD
!!CSD       WRITE(8,*) radius(l),NucleationScaleFactor*nCrystal(l,3,20,1,1)
!!CSD    END DO

!!CSD  CLOSE(UNIT=8,STATUS='keep')

!!CSD  fn='CSD4'
!!CSD  ilength = 4
!!CSD  CALL newfile(fn,suf1,fnv,nint,ilength)
!!CSD  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
!!CSD  WRITE(8,2283) PrintTime
!!CSD  WRITE(8,111)
!!CSD  WRITE(8,*)  ' Radius    nCrystals'

!!CSD    DO l = 1,nCSD
!!CSD       WRITE(8,*) radius(l),nCrystal(l,4,2,1,1)
!!CSD    END DO

!!CSD  CLOSE(UNIT=8,STATUS='keep')

!!CSD  fn='LinearGrowth4'
!!CSD  ilength = 13
!!CSD  CALL newfile(fn,suf1,fnv,nint,ilength)
!!CSD  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
!!CSD  WRITE(8,2283) PrintTime
!!CSD  WRITE(8,111)
!!CSD  WRITE(8,*)  ' Radius    nCrystals'

!!CSD    DO l = 1,nCSD
!!CSD       WRITE(8,*) radius(l),LinearGrowthRate(l,k,10,jy,jz)
!!CSD    END DO

!!CSD  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (nrct > 0) THEN

fn = 'porosity'
ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,112)
112 FORMAT('# Units: % Porosity')
jy = 1
jz = 1
DO jx = 1,nx
  porprt = por(jx,jy,jz)*100.0
  WRITE(8,184) x(jx)*OutputDistanceScale,porprt
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the saturation indices of the minerals (log Q/K).

fn='saturation'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,113)
113 FORMAT('# Units: Dimensionless (Log10(Q/Keq)')
WRITE(8,2285)  (uminprnt(k),k=1,nrct)
jy = 1
jz = 1
DO jx = 1,nx
!!  CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt)
  CALL satcalc(ncomp,nrct,jx,jy,jz)
  DO k = 1,nrct
    dsat(k) = silog(1,k)
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(dsat(k),k=1,nrct)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')
END IF

IF (KateMaher) THEN
  fn='IsotopeRatio'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,2289) 
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidRatio = 18246.5*muUranium234Bulk(jx,jy,jz)/muUranium238Bulk(jx,jy,jz)
    FluidRatio = 18050.5*s(ik234U,jx,jy,jz)/s(ik238U,jx,jy,jz)
	Solid234U = muUranium234Bulk(jx,jy,jz)
	Solid238U = muUranium238Bulk(jx,jy,jz)
	SolidCa = muCalciumBulk(jx,jy,jz)
    SolidRatioMarine1 = 18246.5*mumin(1,kMarineCalcite,ik234U)/mumin(1,kMarineCalcite,ik238U)
    SolidRatioMarine2 = 18246.5*mumin(2,kMarineCalcite,ik234U)/mumin(2,kMarineCalcite,ik238U)
    WRITE(8,184) x(jx)*OutputDistanceScale,SolidRatio,FluidRatio, Solid234U, Solid238U, SolidCa,SolidRatioMarine1,SolidRatioMarine2
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='temperature'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,114)
  114 FORMAT('# Units: Degrees C')
  WRITE(8,2289) 
  jy = 1
  jz = 1
  DO jx = 1,nx
    WRITE(8,184) x(jx)*OutputDistanceScale,t(jx,jy,jz)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF  

!  Write out pressure

IF (calculateflow) THEN
  fn='pressure'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,115)
  115 FORMAT('# Units: Pascals')
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,pres(jx,jy,jz)
    END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!  Write out Darcy fluxes

fn='velocity'
ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,116)
116 FORMAT('# Units: m^3 fluid/m^2 porous medium/yr')
jy = 1
jz = 1
DO jx = 1,nx
  WRITE(8,184) x(jx)*OutputDistanceScale,qx(jx,1,1)
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out pickup file

!  Write out the alkalinity

IF (giambalvo) THEN

fn='alkalinity'
ilength = 10
CALL NEWFILE(fn,suf1,fnv,nint,ilength)
OPEN(unit=8,file=fnv,access='sequential',status='unknown')
WRITE(8,2283) realtime
WRITE(8,2287)
DO jx = 1,nx
  alk = -s(ikph,jx,jy,jz)
  WRITE(8,184) x(jx)*OutputDistanceScale,alk
END DO
CLOSE(unit=8,status='keep')

!  Write out the advective and diffusive fluxes across top boundary (jx=1)
!  and bottom boundary (jx=nx).
!  The calculation (in fx.f) does not work when dispersion .ne. 0)
fn='flux'
ilength = 4
CALL newfile (fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv,ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) realtime
WRITE(8,2296) netflowx(0,1,1)
WRITE(8,2290)'adv top','dif top','tot top','adv bot',  &
    'dif bot','tot bot', 'Tb / Tt', 'nor top', 'nor bot'
WRITE(8,2291)
jy = 1
jz = 1
DO ik = 1, ncomp
  tflux_top = advflux_x(ik,jy,1)+dflux_x(ik,jy,1)
  tflux_bot = advflux_x(ik,jy,2)+dflux_x(ik,jy,2)
  top_norm = tflux_top / ABS(advflux_x(ik,jy,2))
  bot_norm = tflux_bot / ABS(advflux_x(ik,jy,2))
  WRITE(8,2292) ulabprnt(ik), advflux_x(ik,jy,1), dflux_x(ik,jy,1),  &
      tflux_top, advflux_x(ik,jy,2), dflux_x(ik,jy,2),  &
      tflux_bot, tflux_top/tflux_bot, top_norm, bot_norm
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the advective and diffusive fluxes across top boundary (jx=1)
!  with advective relative to seawater concentrations
!  The calculation (in fx.f) does not work when dispersion .ne. 0)
!  Because upflow is negative, >0 is SINK to ocean, <0 is SOURCE to ocean
!  Also write normalized flux (net out sediment/basement direct to ocean)

fn='netflux'
ilength = 7
CALL newfile (fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv,ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) realtime
WRITE(8,2297) netflowx(0,1,1)
WRITE(8,2298) netflowx(nx,1,1)
WRITE(8,2290)'NetAdvT','Dif top','NetTotT', 'NetAdvB','NTT/NAB','ConcSW'
WRITE(8,2293)
jy = 1
jz = 1
DO ik = 1, ncomp
  aflux_net = advflux_x(ik,jy,1) - sbnd(ik,1)*netflowx(0,jy,jz) *ro(1,jy,jz)
  tflux_top = aflux_net + dflux_x(ik,jy,1)
  
!  normalized to advective flux out of basement directly to ocean
  ad_net_bot = advflux_x(ik,jy,2) - sbnd(ik,1) * netflowx(nx,jy,jz)*ro(nx,jy,jz)
  top_norm = tflux_top / ad_net_bot
    WRITE(8,2294) ulabprnt(ik), aflux_net, dflux_x(ik,jy,1),  &
        tflux_top, ad_net_bot, top_norm, sbnd(ik,1)
END DO
CLOSE(UNIT=8,STATUS='keep')

END IF

502 FORMAT('temperature    ' ,f8.2)
503 FORMAT(a20,4X,1PE12.4)
504 FORMAT('END')
182 FORMAT(80(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE16.8))

!2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2283 FORMAT('# Time      ',2X,1PE12.4)
2284 FORMAT('    Distance ',a18)
2282 FORMAT('    Distance ','        pH')
2281 FORMAT('    Distance ',4X,a18)
2285 FORMAT('    Distance       ',100(1X,a16))
2299 FORMAT('    Distance    ',100(1X,a16))
2288 FORMAT('    Distance    ','  pH           ',100(1X,a16))
2286 FORMAT('    Distance    ',100(1X,a16))
2289 FORMAT('    Distance    ',1x,' CalciteAR',1x,' PoreWaterAR', 1x, '234U', 1x, '238U', 1x, 'Ca')
514 FORMAT(1X,a10,1X,1PE11.4)
513 FORMAT(1X,a18,1X,1PE12.4,5X,a12)

600 FORMAT(2X,f10.2,2X,a15)
201 FORMAT(2X,a18,2X,f8.2)
202 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,1PE12.3,2x,a8)
211 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,'            ',2x,a8)
203 FORMAT(2X,a18)
207 FORMAT('                 ','          Log',1X,'       Log',1x,  &
'              ',1x,'          ', '       Activity')
204 FORMAT(' Species         ','     Molality',1X,'  Activity',1x,  &
'     Molality ',1x,'    Activity ', '  Coefficient','    Type')
205  FORMAT(2X,'Total Charge    = ',1pe10.3)

509 FORMAT(2X,a18,2X,f12.4)
510 FORMAT(2X,'GEOCHEMICAL CONDITION NUMBER',i3)
555 FORMAT(2X,'Temperature (C) = ',f10.3)
411 FORMAT(2X,'Ionic Strength  = ',f10.3)
5022 FORMAT(2X,'Solution pH     = ',f10.3)
5023 FORMAT(2X,'Solution pe     = ',f10.3)
5024 FORMAT(2X,'Solution Eh     = ',f10.3)
512 FORMAT(' Basis species    ','     Molality  ', '  Constraint type')
511 FORMAT(1X,a18,1X,1PE12.4,5X,a13,1X,a18)

542  FORMAT(5X,'Primary species ',a18,1X,'at grid pt. ',i5)
544  FORMAT(5X,'Check geochem. condition # ',i2)
543  FORMAT(5X,'Constraint mineral ',2X,a18)
643  FORMAT(5X,'Constraint gas',2X,a18)
2287 FORMAT('#   Distance ',' alkalinity (eq/kg)')
2290 FORMAT('#  Component',5X,12(a7,7X))
2291 FORMAT('#           ',5X,6('mol/m2/yr',5X))
2292 FORMAT(a15, 9(1X,1PE13.6), 2(1X,1I2), 2(1X,1PE13.6))
2293 FORMAT('#           ',5X,3('mol/m2/yr',5X))
2294 FORMAT(a15, 9(1X,1PE13.6))
2296 FORMAT('# Net flow at top: ',1X,1PE13.6)
2297 FORMAT('# Net flow at top:    ',1X,1PE13.6)
2298 FORMAT('# Net flow at bottom: ',1X,1PE13.6)



RETURN
END SUBROUTINE GraphicsKaleidagraph
!  *******************************************************
