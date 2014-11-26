SUBROUTINE reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium
USE temperature
USE runtime
USE NanoCrystal
USE isotope

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(IN)                                        :: nkin
INTEGER(I4B), INTENT(IN)                                        :: nrct
INTEGER(I4B), INTENT(IN)                                        :: nspec
INTEGER(I4B), INTENT(IN)                                        :: nexchange
INTEGER(I4B), INTENT(IN)                                        :: nsurf
INTEGER(I4B), INTENT(IN)                                        :: ndecay
INTEGER(I4B), INTENT(IN)                                        :: jx
INTEGER(I4B), INTENT(IN)                                        :: jy
INTEGER(I4B), INTENT(IN)                                        :: jz

REAL(DP), INTENT(IN)                                            :: delt
REAL(DP), INTENT(IN)                                            :: time

!  Internal variables and arrays

REAL(DP)                                                        :: tk
REAL(DP)                                                        :: tkinv
REAL(DP)                                                        :: reft
REAL(DP)                                                        :: porfactor
REAL(DP)                                                        :: sumiap
REAL(DP)                                                        :: term2
REAL(DP)                                                        :: term1
REAL(DP)                                                        :: vcheck 
REAL(DP)                                                        :: check 
REAL(DP)                                                        :: checkmonod
REAL(DP)                                                        :: sign

INTEGER(I4B)                                                    :: k
INTEGER(I4B)                                                    :: np
INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: is
INTEGER(I4B)                                                    :: kk
INTEGER(I4B)                                                    :: nex
INTEGER(I4B)                                                    :: ixx
INTEGER(I4B)                                                    :: l
INTEGER(I4B)                                                    :: kJapaneseGlassNoSilica

REAL(DP)                                                        :: RateFactor
REAL(DP)                                                        :: CubicTerm
REAL(DP)                                                        :: CorrectedTerm
!!REAL(DP)                                                        :: Kformation
REAL(DP)                                                        :: Al

REAL(DP)                                                        :: dependTMP
REAL(DP)                                                        :: surftmp
REAL(DP)                                                        :: MinConvert
REAL(DP)                                                        :: CheckRgas
REAL(DP)                                                        :: AffinityTerm

REAL(DP)                                                        :: Astar
REAL(DP)                                                        :: Bstar
REAL(DP)                                                        :: Sstar
REAL(DP)                                                        :: denominator

!!!REAL(DP)                                                        :: MoleFraction40
!!!REAL(DP)                                                        :: MoleFraction44
!!!REAL(DP)                                                        :: MoleFraction32
!!!REAL(DP)                                                        :: MoleFraction34
!!!REAL(DP)                                                        :: MoleFraction32S
!!!REAL(DP)                                                        :: MoleFraction34S

!!!REAL(DP)                                                        :: RateFraction
!!!REAL(DP)                                                        :: MoleFraction44Mineral
!!!REAL(DP)                                                        :: MoleFraction40Mineral
!!!REAL(DP)                                                        :: CurrentMoleFraction40
!!!REAL(DP)                                                        :: CurrentMoleFraction44
!!!REAL(DP)                                                        :: CheckMoleFraction
!!!REAL(DP)                                                        :: CheckLogQ
!!!REAL(DP)                                                        :: CheckLogQ10
!!!REAL(DP)                                                        :: kBackward
!!!REAL(DP)                                                        :: kForward
!!!REAL(DP)                                                        :: rminCheck

REAL(DP)                                                        :: MoleFractionMineral

INTEGER(I4B)                                                    :: kIsotopologue
INTEGER(I4B)                                                    :: Isotopologue
INTEGER(I4B)                                                    :: kMineralRare
INTEGER(I4B)                                                    :: kMineralCommon
INTEGER(I4B)                                                    :: iPrimaryRare
INTEGER(I4B)                                                    :: iPrimaryCommon

REAL(DP)                                                        :: SaturationCloseToEquilibrium
REAL(DP)                                                        :: Term1RateConstant
REAL(DP)                                                        :: CheckRateAt1
REAL(DP)                                                        :: CheckRateAt2
REAL(DP)                                                        :: CallItEquilibrium

REAL(DP)                                                        :: CalciumCarbonateRatioEffect

! biomass
INTEGER(I4B)                                                    :: ib, jj
! end biomass

!!LOGICAL(LGT)                                                    :: UseDissolutionOnly
!!LOGICAL(LGT)                                                    :: checkSaturation

!! CrystalSizeDistribution
REAL(DP)                                                        :: SizeDepKeq
REAL(DP)                                                        :: AdjustKeq

rmin = 0.0d0

!!CSD    sigma(3) = 200.0/1000.0
!!CSD    sigma(4) = 90.0/1000.0

!  This routine calculates the reaction rates of the individual
!  minerals.

tk = t(jx,jy,jz) + 273.15D0
tkinv = 1.0D0/tk
reft = 1.0D0/298.15D0

porfactor = (por(jx,jy,jz)/porin(jx,jy,jz))**(0.666666D0)
porfactor = 1.0d0
!!porfactor = 1.0d0
!!porfactor = por(jx,jy,jz)
IF (DampRateInLowPorosity .AND. por(jx,jy,jz) < 0.001) THEN
  porfactor = porfactor*PorosityDamp
END IF

IF (nIsotopePrimary > 0) THEN
 UseAqueousMoleFraction = .FALSE.
END IF

!!porfactor = 1.0d0

!! For isotopes, calculate mole fractions based on aqueous geochemistry (for No Back Reaction case)

DO kIsotopologue = 1,nIsotopeMineral
  kMineralRare = kIsotopeRare(kIsotopologue)
  KMineralCommon = kIsotopeCommon(kIsotopologue)
  isotopologue = PointerToPrimaryIsotope(kIsotopologue)
  iPrimaryCommon = isotopeCommon(Isotopologue)
  IF (IsotopeBackReactionOption(kIsotopologue) == 'bulk' .OR. time < LagSurface/365.0) THEN
    denominator = volfx(kMineralRare,jx,jy,jz) + volfx(kMineralCommon,jx,jy,jz) 
    IF (denominator == 0.0d0 .OR. time < LagSurface/365.0) THEN
      UseAqueousMoleFraction(kIsotopologue) = .TRUE.
    ELSE
      MoleFractionMineralRare(kIsotopologue) = (volfx(kMineralRare,jx,jy,jz)/denominator)
      MoleFractionMineralCommon(kIsotopologue) = (volfx(kMineralCommon,jx,jy,jz)/denominator)
    END IF
  ELSE IF (IsotopeBackReactionOption(kIsotopologue) == 'surface') THEN
    denominator = volsave(kMineralRare,jx,jy,jz) + volsave(kMineralCommon,jx,jy,jz)
    IF (denominator == 0.0d0 .OR. time<LagSurface/365.0) THEN
      UseAqueousMoleFraction(kIsotopologue) = .TRUE.
    ELSE
      MoleFractionMineralRare(kIsotopologue) = volsave(kMineralRare,jx,jy,jz)/denominator
      MoleFractionMineralCommon(kIsotopologue) = volsave(kMineralCommon,jx,jy,jz)/denominator
    END IF
  ELSE
    UseAqueousMoleFraction(kIsotopologue) = .TRUE.
  END IF
END DO

DO Isotopologue = 1,nIsotopePrimary
  iPrimaryRare = isotopeRare(Isotopologue)
  iPrimaryCommon = isotopeCommon(Isotopologue)
  denominator = sp10(iPrimaryRare,jx,jy,jz) + sp10(iPrimaryCommon,jx,jy,jz)
  MoleFractionAqueousRare(Isotopologue) = (sp10(iPrimaryRare,jx,jy,jz)/denominator)
  MoleFractionAqueousCommon(Isotopologue) = (sp10(iPrimaryCommon,jx,jy,jz)/denominator)
END DO


!! For dissolution, use the mole fraction of 44Ca and 40Ca in the mineral from previous time step

!!!!!!  IF (volfx(6,jx,jy,jz) == 0.0d0 .AND. volfx(5,jx,jy,jz) == 0.0d0 .OR. UseAqueousMoleFraction) THEN

!!!!!!    IF (time < LagSurface/365.0 .OR. UseBulkMineral) THEN
!!!!!      MoleFraction44Mineral = volfx(6,jx,jy,jz)/( volfx(6,jx,jy,jz) + volfx(5,jx,jy,jz) )

!!!!!      MoleFraction44Mineral = VolSave(6,jx,jy,jz)/( VolSave(6,jx,jy,jz) + VolSave(5,jx,jy,jz) )


decay_correct = 1.0D0
snorm = 0.0d0
ivolume = 0
checkSaturation = .FALSE.
UseDissolutionOnly = .FALSE.

DO k = 1,nkin

  checkSaturation = .FALSE.
  UseDissolutionOnly = .FALSE.

  dppt(k,jx,jy,jz) = 0.0D0

!!  **********  CrystalSizeDistribution  *********************************************
!!CSD  IF (CrystalSizeDistribution(k)) THEN  
!!CSD    DO l = 1,nCSD
!!CSD      AdjustKeq = 2.0d0*sigma(k)*volmol(k)/(rgas*Tk)* (1.0d0/radius(l))
!!CSD      AdjustKeqLog(l,k) = DLOG(AdjustKeq)
!!CSD    END DO
!!CSD  END IF
!!  **********  End of CrystalSizeDistribution  **************************************  

  DO np = 1,nreactmin(k)
    
    IF (t(jx,jy,jz) == 25.0D0) THEN
      actenergy(np,k) = 1.0D0
    ELSE
      actenergy(np,k) = DEXP( (ea(np,k)/rgasKCAL)*(reft-tkinv) )
      CheckRgas = rgasKCAL
    END IF
    
    IF (KateMaher) THEN
      mumin(np,kUCalcite,ikCa) = muCalcium(jx,jy,jz)
      mumin(np,kUCalcite,ik234U) = muUranium234(jx,jy,jz)
      mumin(np,kUCalcite,ik238U) = muUranium238(jx,jy,jz)
    END IF

!!  Imintype = 1 --> TST
!!  Imintype = 2 --> Monod kinetics (NOTE: Should add the affinity term here)
!!  Imintype = 3 --> Irreversible kinetics
!!  Imintype = 4 --> Precipitation only, thermodynamic term
!!  Imintype = 5 --> Dissolution only, thermodynamic term
!!  Imintype = 6 --> Forward one-way, surface area = 1.0 while mineral is present
!!  Imintype = 7 --> Reverse one-way, surface area = 1.0 while mineral is present
!!  Imintype = 8 --> Microbially mediated monod with theormodynamic factor, and biomass
!!  Imintype = 9 --> biomass first order decay !!smr-2011-04-29

!!  Calculate thermodynamic term
    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN  !! Simple Monod (Option 2 with no Ft) or Irreversible (Option 3)

      silog(np,k) = -500.0D0
      si(np,k) = 0.0D0

    ELSE

!!  Cross-Affinity Case
      IF (kcrossaff(np,k) /= 0) THEN

        kk = kcrossaff(np,k)

        IF (ikh2o /= 0) THEN

          sumiap = 0.0D0
          DO i = 1,ncomp
            IF (ulab(i) == 'H2O') THEN
              sumiap = sumiap + decay_correct(i,k)*mumin(1,kk,i)*(gam(i,jx,jy,jz))
            ELSE
              sumiap = sumiap + decay_correct(i,k)*mumin(1,kk,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
            END IF
          END DO
  
        ELSE

          sumiap = 0.0D0
          DO i = 1,ncomp
            sumiap = sumiap + decay_correct(i,k)*mumin(1,kk,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
          END DO

        END IF

        silog(np,k) = (sumiap - keqmin(1,kk,jx,jy,jz))/clg
        si(np,k) = 10**(silog(np,k))

!!  Biomass case
      ELSE IF (imintype(np,k) == 8) THEN             !! Microbially mediated, with thermodynamic factor Ft

        jj = p_cat_min(k)

        IF (ikh2o /= 0) THEN

          sumiap = 0.0D0
          DO i = 1,ncomp
            IF (ulab(i) == 'H2O') THEN
              sumiap = sumiap + muminTMP(np,jj,i)*(gam(i,jx,jy,jz))
            ELSE
              sumiap = sumiap + muminTMP(np,jj,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
            END IF
          END DO

        ELSE

          sumiap = 0.0D0
          DO i = 1,ncomp
            sumiap = sumiap + muminTMP(np,jj,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
          END DO

        END IF

        silog(np,k) = (sumiap - keqminTMP(np,jj) - BQ_min(np,jj)/(rgas*Tk))/clg    !!  BQ in kJ/e-mole
        siln(np,k) = clg*silog(np,k)
        si(np,k) = 10**(silog(np,k))

!!  Base case
      ELSE    
            
        IF (ikh2o /= 0) THEN     
                     
          sumiap = 0.0D0
          DO i = 1,ncomp
            IF (ulab(i) == 'H2O') THEN
              sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(gam(i,jx,jy,jz))
            ELSE
              sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
            END IF
          END DO

        ELSE

          sumiap = 0.0D0
          DO i = 1,ncomp
            sumiap = sumiap + decay_correct(i,k)*mumin(np,k,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
          END DO

        END IF
   
        IF (nIsotopeMineral > 0) THEN

!!  Now modify the Ion Activity Product (Q) for the case of a solid solution (activity /= 1)

        IF (IsotopeMineralRare(k)) THEN

          kIsotopologue = kPointerIsotope(k)

          kMineralRare = kIsotopeRare(kIsotopologue)
          KMineralCommon = kIsotopeCommon(kIsotopologue)
          isotopologue = PointerToPrimaryIsotope(kIsotopologue)
          iPrimaryCommon = isotopeCommon(Isotopologue)

          IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN

              MoleFractionMineral = MoleFractionAqueousRare(isotopologue)

          ELSE

            IF (MoleFractionMineralRare(kIsotopologue) == 0.0d0) THEN
              MoleFractionMineral = 1.0d0
            ELSE

                MoleFractionMineral = MoleFractionMineralRare(isotopologue)


            END IF
          END IF

          sumiap = sumiap - (mumin(1,kMineralCommon,iPrimaryCommon))*DLOG(MoleFractionMineral)

        ELSE IF (IsotopeMineralCommon(k)) THEN

          kIsotopologue = kPointerIsotope(k)

          kMineralRare = kIsotopeRare(kIsotopologue)
          KMineralCommon = kIsotopeCommon(kIsotopologue)
          isotopologue = PointerToPrimaryIsotope(kIsotopologue)
          iPrimaryCommon = isotopeCommon(Isotopologue)

          IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN

                MoleFractionMineral = MoleFractionAqueousCommon(isotopologue)

          ELSE
            IF (MoleFractionMineralCommon(kIsotopologue) == 0.0d0) THEN
              MoleFractionMineral = 1.0d0
            ELSE

                MoleFractionMineral = MoleFractionMineralCommon(isotopologue)

            END IF
          END IF

          sumiap = sumiap - (mumin(1,kMineralCommon,iPrimaryCommon))*DLOG(MoleFractionMineral)

        ELSE
          CONTINUE
        END IF

        END IF   !! Block where nIsotopeMineral > 0

        silog(np,k) = (sumiap - keqmin(np,k,jx,jy,jz))/clg
        siln(np,k) = clg*silog(np,k)
        si(np,k) = 10**(silog(np,k))
        
!!CSD!!  **********  CrystalSizeDistribution  **********************************
!!CSD        IF (CrystalSizeDistribution(k)) THEN
!!CSD
!!CSD!!        Cycle over crystal size distribution
!!CSD
!!CSD          DO l = 1,nCSD
!!CSD            SizeDepKeq = keqmin(np,k,jx,jy,jz) + AdjustKeqLog(l,k)
!!CSD            silogCrystalSize(l,np,k) =  (sumiap - SizeDepKeq)/clg
!!CSD!!          Assumes half sphere--> 0.5 * 4*pi*r^2
!!CSD            AreaCrystalSize(l,k) = NucleationScaleFactor*DFLOAT(nCrystal(l,k,jx,jy,jz)) * 2.0d0*pi*radius(l)*radius(l)
!!CSD            siCrystalSize(l,np,k) = 10**(silogCrystalSize(l,np,k)) 
!!CSD          END DO
!!CSD   
!!CSD        END IF
!!CSD!!  **********  CrystalSizeDistribution  *********************************************

    
!!        IF (DABS(silog(np,k)) < 0.00001 .AND. checkSaturation) THEN
!!          UseDissolutionOnly = .TRUE.
!!        END IF

      END IF

    END IF

!!  Decide on the value of the reactive surface area depending on
!!  whether the mineral is under or super-saturated and whether
!!  it has a zero volume fraction.

!!  Reset the surface areas for minerals that are associated with other minerals
!!  "MineralID" is the pointer to the mineral number whose volume fraction is being tracked
    
    IF (silog(np,k) >= 0.0D0 .AND. iarea(k,jinit(jx,jy,jz)) == 0) THEN       !! Supersaturated AND bulk_surface_area option

!!  Associate mineral with another mineral (surface area and volume fraction)
      IF (MineralAssociate(k)) THEN

        IF (MineralID(k) < k) THEN                  !!  NOTE: This requires that the mineral that is associated with is earlier in list
          surf(k) = surf(MineralID(k))
        ELSE
          surf(k) = areain(MineralID(k),jinit(jx,jy,jz))*porfactor
        END IF

      ELSE

!! NOTE:  Perhaps should change this to a specific option for supersaturation, with default = .TRUE.  ??
        IF (SetSurfaceAreaConstant) THEN
          surf(k) = areain(k,jinit(jx,jy,jz))*porfactor
        ELSE
          surf(k) = area(k,jx,jy,jz)*porfactor            
        END IF

        surf(k) = areain(k,jinit(jx,jy,jz))*porfactor

      END IF

!!  Case where either undersaturated OR using the specific_surface_area option  

    ELSE                                                

      IF (MineralAssociate(k)) THEN

        if (MineralID(k) < k) then
          surf(k) = surf(MineralID(k))
        else

          IF (porfactor < 0.01d0) THEN
            surf(k) = area(MineralID(k),jx,jy,jz)*porfactor
          ELSE
            surf(k) = area(MineralID(k),jx,jy,jz)
          END IF

        end if

      ELSE

        IF (porfactor < 0.01d0) THEN
          surf(k) = area(k,jx,jy,jz)*porfactor
        ELSE
          surf(k) = area(k,jx,jy,jz)
        END IF

!!!        IF (SetSurfaceAreaConstant) THEN
!!!          surf(k) = areain(k,jinit(jx,jy,jz))*porfactor
!!!        ELSE
!!!          surf(k) = area(k,jx,jy,jz)*porfactor
!!!        END IF

      END IF

    END IF

!! Reset surface area for irreversible or Monod reaction
    IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN

      IF (MineralAssociate(k)) THEN
        surf(k) = area(MineralID(k),jx,jy,jz)
      ELSE
        surf(k) = area(k,jx,jy,jz)
      END IF

    END IF

    IF (KateMaher) THEN
      IF (k==KUcalcite .AND. rlabel(np,k) == 'recrystallize') THEN
        surf(k) = area(kMarineCalcite,jx,jy,jz)
      END IF
    END IF
    
    IF (surf(k) == 0.0D0) THEN

      rmin(np,k) = 0.0D0
      dppt(k,jx,jy,jz) = 0.0D0
      pre_rmin(np,k) = 1.0D0

    ELSE
          
!!    TST, irreversible, or PrecipitationOnly
      IF (imintype(np,k) == 1 .OR. imintype(np,k) == 3 .OR. imintype(np,k) == 4 .OR. imintype(np,k) == 5) THEN 
        
        term2 = 0.0D0
        DO kk = 1,ndepend(np,k)
          i = idepend(kk,np,k)
          IF (depend(kk,np,k) < 0.000) THEN
            IF (itot_min(kk,np,k) == 1) THEN
              term2 = term2 + depend(kk,np,k)*DLOG(s(i,jx,jy,jz))
            ELSE
              IF (ulab(i) == 'H2O') THEN
                term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz))
              ELSE
                term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz)+sp(i,jx,jy,jz))
              END IF
            END IF
          ELSE
            IF (itot_min(kk,np,k) == 1) THEN
              term2 = term2 + depend(kk,np,k)*DLOG(s(i,jx,jy,jz))
            ELSE
              IF (ulab(i) == 'H2O') THEN
                term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz))
              ELSE
                term2 = term2 + depend(kk,np,k)*(gam(i,jx,jy,jz)+sp(i,jx,jy,jz))
              END IF
            END IF
          END IF
        END DO
        DO kk = 1,ndependex(np,k)
          nex = ixdepend(kk,np,k)
          ixx = ixlink(nex)
          term2 = term2 + dependex(kk,np,k)*spex(nex+nexchange,jx,jy,jz)
        END DO
        DO kk = 1,ndependsurf(np,k)
          is = isdepend(kk,np,k)
          dependTMP = dependsurf(kk,np,k)
          IF (is <= nsurf) THEN
            term2 = term2 + dependsurf(kk,np,k)*spsurf(is,jx,jy,jz)
          ELSE
            term2 = term2 + dependTMP*DLOG(spsurf10(is,jx,jy,jz))        
          END IF
        END DO

        IF (HyperbolicInhibition(np,k)) THEN
!!          Kformation = 0.0000149823D0
          Al = sp10(HyperbolicInhibitionPointer(np,k),jx,jy,jz)
          denominator = Kformation(np,k) + Al**HyperbolicInhibitionDepend(np,k)
          term2 = Kformation(np,k)/denominator
        END IF
        
!!        ***** * Monod kinetics (without thermodynamic term)  **********
        
      ELSE IF (imintype(np,k) == 2) THEN     !  Monod rate law

!!  If in Monod rate law there is a Monod term for the mineral itself, do not use the mineral surface area (reset to 1.0) and use the hyperbolic expression  

        term2 = 1.0D0
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            surf(k) = 1.0d0                                                        !! Resets surface area to a constant 1
            MinConvert = volfx(k,jx,jy,jz)/(volmol(k)*por(jx,jy,jz)*ro(jx,jy,jz))  !! Converts mineral volume fraction to moles mineral per kg fluid (molality)                                  
            checkmonod =  MinConvert/(MinConvert+halfsat(kk,np,k))
            term2 = term2 * checkmonod 
          ELSE                                                                     !! Monod term depends on aqueous species
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN
              checkmonod =  s(i,jx,jy,jz)/(s(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod             
            ELSE
              checkmonod =  sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod
            END IF
          END IF
        END DO
        DO kk = 1,ninhibit(np,k)

          i = inhibit(kk,np,k)
          IF (itot_inhibit(kk,np,k) == 1) THEN
            term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+s(i,jx,jy,jz))
          ELSE
            term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+sp10(i,jx,jy,jz))
          END IF
        END DO
!!        IF (umin(k)=='JapaneseGlassSilicaOnly') THEN 
!!          kJapaneseGlassNoSilica = 2                                                     
!!          checkmonod =  0.01d0/(volfx(kJapaneseGlassNoSilica,jx,jy,jz) + 0.01d0)
!!          term2 = term2 * checkmonod 
!!        END IF
      ELSE IF (imintype(np,k) == 6 .OR. imintype(np,k) == 7) THEN  !! Forward or reverse, with surface area set constant at 1.0

        CONTINUE

!!    Monod rate law with biomass and thermodynamic factor, F_T     
      ELSE IF (imintype(np,k) == 8) THEN     

        jj = p_cat_min(k)

!!      Metabolic lag function (see Wood et al, 1995)
        IF (UseMetabolicLagMineral(np,jj)) THEN

!!        S* == Critical substrate concentration
!!        A* == Metabolic lag (in years)
!!        B* == Metabolic lag + ramp up period (so B* - A* = ramp up period (yrs)
          Astar = LagTimeMineral(np,jj)/365.0d0         !! Days converted to years
          Bstar = RampTimeMineral(np,jj)/365.0d0        !! Days converted to years
          Sstar = ThresholdConcentrationMineral(np,jj)  !! Critical concentration of substrate (acetate in this case)

          IF (sn(SubstrateForLagMineral(np,jj),jx,jy,jz) > Sstar .AND. tauZeroMineral(np,jj,jx,jy,jz) == 0.0d0) THEN
            tauZeroMineral(np,jj,jx,jy,jz) = time
          END IF

          IF (tauZeroMineral(np,jj,jx,jy,jz) == 0.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0d0
          ELSE IF (tauZeroMineral(np,jj,jx,jy,jz) > 0.0d0 .AND. time < (tauZeroMineral(np,jj,jx,jy,jz)+Astar) ) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0D0
          ELSE    
            denominator = (tauZeroMineral(np,jj,jx,jy,jz)+Astar+Bstar) - (tauZeroMineral(np,jj,jx,jy,jz)+Astar)
            IF (denominator == 0.0d0) THEN
              MetabolicLagMineral(np,jj,jx,jy,jz) = 1.0d0
            ELSE
              MetabolicLagMineral(np,jj,jx,jy,jz) =  (time -(tauZeroMineral(np,jj,jx,jy,jz)+Astar) )/denominator
            END IF
          END IF
          IF (MetabolicLagMineral(np,jj,jx,jy,jz) >= 1.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 1.0D0
          ELSE IF (MetabolicLagMineral(np,jj,jx,jy,jz) <= 0.0d0) THEN
            MetabolicLagMineral(np,jj,jx,jy,jz) = 0.0D0
          ELSE
            CONTINUE
          END IF

        END IF

!!  If in Monod rate law there is a Monod term for the mineral itself, do not use the mineral surface area (reset to 1.0) and use the hyperbolic expression  
             
        term2 = 1.0D0
        DO kk = 1,nmonod(np,k)
          IF (imonod(kk,np,k) == 0 .AND. kmonod(kk,np,k) == 1) THEN                !! Mineral Monod reaction depends on its own concentration
            surf(k) = 1.0d0                                                        !! Resets surface area to a constant 1
            MinConvert = volfx(k,jx,jy,jz)/(volmol(k)*por(jx,jy,jz)*ro(jx,jy,jz))  !! Converts mineral volume fraction to moles mineral per kg fluid (molality)                                  
            checkmonod =  MinConvert/(MinConvert+halfsat(kk,np,k))
            term2 = term2 * checkmonod 
          ELSE                                                                     !! Monod term depends on aqueous species
            i = imonod(kk,np,k)
            IF (itot_monod(kk,np,k) == 1) THEN
              checkmonod =  s(i,jx,jy,jz)/(s(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod             
            ELSE
              checkmonod =  sp10(i,jx,jy,jz)/(sp10(i,jx,jy,jz)+halfsat(kk,np,k))
              term2 = term2 * checkmonod
            END IF
          END IF
        END DO

! frankfurt - sergi - add inhibition
        !DO kk = 1,ninhibit(np,k)
        !  i = inhibit(kk,np,k)
        !  IF (itot_inhibit(kk,np,k) == 1) THEN
        !    term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+s(i,jx,jy,jz))
        !  ELSE
        !    term2 = term2 * rinhibit(kk,np,k)/(rinhibit(kk,np,k)+sp10(i,jx,jy,jz))
        !  END IF
        ! END DO        
! biomass end

! sergi: first order biomass decay - may 2011
      
      ELSE IF (imintype(np,k) == 9) THEN

        term2 = volfx(biomass_decay(np,nkin),jx,jy,jz)
        !!write(*,*)'sergi: volume of biomass: ',volfx(biomass_decay(np,nkin),jx,jy,jz)
        pre_rmin(np,k) = term2
        continue
        
! end sergi: first order biomass decay - may 2011

      ELSE

        WRITE(*,*)
        WRITE(*,*) ' Rate formulation not recognized'
        WRITE(*,*) ' Imintype = ', imintype(np,k)
        WRITE(*,*) ' For mineral = ',umin(k)
        WRITE(*,*) ' Parallel rate law = ', np
        WRITE(*,*)
        READ(*,*)
        STOP

      END IF
      
      IF (imintype(np,k) == 2) THEN     ! Monod kinetics

        pre_rmin(np,k) = term2

! biomass
      ELSE IF (imintype(np,k) == 8) THEN   ! MonodBiomass kinetics

        pre_rmin(np,k) = term2
! biomass end        

! sergi: first order biomass decay - may 2011
      ELSE IF (imintype(np,k) == 9) THEN 

        pre_rmin(np,k) = term2
! end sergi: first order biomass decay - may 2011

      ELSE IF (HyperbolicInhibition(np,k)) THEN  
        pre_rmin(np,k) = term2
      ELSE

        IF (term2 == 0.0D0) THEN
          pre_rmin(np,k) = 1.0D0
        ELSE
          pre_rmin(np,k) = DEXP(term2)
        END IF

      END IF
      
!***************** Dependence on mineral saturation state *************
     
      IF (si(np,k) > 1.0D0) THEN
        sign = 1.0D0
      ELSE
        sign = -1.0D0
      END IF

      IF (KateMaher) THEN
        IF (k == KUCalcite .AND. rlabel(np,k) == 'recrystallize') THEN
          sign = 1.0D0
        END IF
      END IF

      IF (imintype(np,k) == 3 .OR. imintype(np,k) == 2) THEN     !!  Monod or irreversible

        snorm(np,k) = 0.0D0

! biomass
      ELSE IF (imintype(np,k) == 8) THEN    !!  Monod, but with thermodynamic factor (Ft) a la Jin and Bethke

! added
        jj = p_cat_min(k)

!!  Add in here the thermodynamic factor, F_T

!!  Formulation assumes that Delta G is negative and that it goes forward, but we have written the ppt reaction with 
!!  the mineral on the left.  So invert the Q/Keq if supersaturated (sign will capture the direction)

        IF ( si(np,k) > 1.0d0) then

          if (direction_min(np,jj) == 1) THEN
            snorm(np,k) = (1.0d0/si(np,k))**(1.0d0/chi_min(np,jj))
          ELSE
            snorm(np,k) = 1.0d0
          END IF

        ELSE
 
          IF( direction_min(np,jj) == -1) THEN
            snorm(np,k) = si(np,k)**(1.0d0/chi_min(np,jj))
          ELSE
            snorm(np,k) = 1.0d0
          END IF

        END IF

        term1 = sign*DABS(snorm(np,k) - 1.0D0)

        IF (direction_min(np,jj) == 1) THEN               !! Precipitation only
          AffinityTerm = MAX(0.0d0,term1)
        ELSE                                         !! Dissolution only
          AffinityTerm = MIN(0.0d0,term1)
        ENDIF
              
! end Monodbiomass

! sergi: first order biomass decay - may 2011
      
      ELSE IF (imintype(np,k) == 9) THEN
      
        term1 = 1.d0
        AffinityTerm = 1.0d0
      
! end sergi: first order biomass decay - may 2011

      ELSE

        IF (AffinityDepend2(np,k) == 1.0D0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)
        ELSE IF (AffinityDepend2(np,k) /= 1.0d0 .AND. AffinityDepend3(np,k) == 1.0d0) THEN
          snorm(np,k) = si(np,k)**AffinityDepend2(np,k)
        ELSE
          snorm(np,k) = DEXP(-AffinityDepend2(np,k)*(DABS(siln(np,k)))**AffinityDepend3(np,k) )
        END IF

      END IF    !! Sergi:  Need to move this down?

      IF (AffinityDepend1(np,k) == 1.0D0) THEN
        term1 = sign*DABS(snorm(np,k) - 1.0D0)
      ELSE      
        term1 = sign*DABS(snorm(np,k) - 1.0D0)**(AffinityDepend1(np,k))
      END IF

      IF (imintype(np,k) == 5 .OR. UseDissolutionOnly) THEN                          !! Dissolution only
!!        IF (UseDissolutionOnly) THEN
!!          AffinityTerm = term1
!!        ELSE
          AffinityTerm = MIN(0.0d0,term1)
!!        END IF
      ELSE IF (imintype(np,k) == 4) THEN                               !! Precipitation only
        AffinityTerm = MAX(0.0d0,term1)
      ELSE
!! General case, which uses (or not) the exponents in "snorm"
        AffinityTerm = term1
      END IF

!!    Sergi:  Need to move the END IF above here??


!************** End of dependence on mineral saturation state *************

! commented out        rmin(np,k) = surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm 

!!    Biomass, with thermodynamic term, Ft         
      IF (imintype(np,k) == 8) THEN   !! MonodBiomass

!       pointer to biomass for current reaction
        jj = p_cat_min(k)
!       pointer to biomass for current reaction
        ib = ibiomass_min(jj)

        IF (UseMetabolicLagMineral(np,jj)) THEN
          rmin(np,k) = MetabolicLagMineral(np,jj,jx,jy,jz)*surf(k)*rate0(np,k)*           &
              actenergy(np,k)*pre_rmin(np,k)*volfx(ib,jx,jy,jz)*AffinityTerm
        ELSE
          rmin(np,k) = surf(k)*rate0(np,k)*           &
             actenergy(np,k)*pre_rmin(np,k)*volfx(ib,jx,jy,jz)*AffinityTerm
        END IF     
! biomass end

      ELSE

        IF (nIsotopeMineral > 0) THEN


          IF (IsotopeMineralCommon(k)) THEN

            kIsotopologue = kPointerIsotope(k)
            kMineralRare = kIsotopeRare(kIsotopologue)
            KMineralCommon = kIsotopeCommon(kIsotopologue)
            isotopologue = PointerToPrimaryIsotope(kIsotopologue)
            iPrimaryCommon = isotopeCommon(Isotopologue)

            IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN
              isotopologue = PointerToPrimaryIsotope(kIsotopologue)
              MoleFractionMineral = MoleFractionAqueousCommon(isotopologue)
            ELSE
              MoleFractionMineral = MoleFractionMineralCommon(kPointerIsotope(k))
            END IF

          ELSE IF (IsotopeMineralRare(k)) THEN

            kIsotopologue = kPointerIsotope(k)
            kMineralRare = kIsotopeRare(kIsotopologue)
            KMineralCommon = kIsotopeCommon(kIsotopologue)
            isotopologue = PointerToPrimaryIsotope(kIsotopologue)
            iPrimaryCommon = isotopeCommon(Isotopologue)

            IF (isotopeBackReactionOption(kIsotopologue) == 'none' .OR. UseAqueousMoleFraction(kIsotopologue)) THEN
              isotopologue = PointerToPrimaryIsotope(kIsotopologue)
              MoleFractionMineral = MoleFractionAqueousRare(isotopologue)

!!            Use the bulk surface area (common)
              surf(k) = surf(kIsotopeCommon(kIsotopologue))
            ELSE
              MoleFractionMineral = MoleFractionMineralRare(kPointerIsotope(k))
!!            Use the bulk surface area (common)
              surf(k) = surf(kIsotopeCommon(kIsotopologue))
            END IF

          ELSE
            MoleFractionMineral = 1.0d0
          END IF
        
        ELSE
          MoleFractionMineral = 1.0d0
        END IF

        rmin(np,k) = MoleFractionMineral*surf(k)*rate0(np,k)*actenergy(np,k)*pre_rmin(np,k)*AffinityTerm
              
      END IF

    END IF

    dppt(k,jx,jy,jz) = dppt(k,jx,jy,jz) + rmin(np,k)
  
  END DO   !  End of npth parallel reaction
  
  IF (MineralAssociate(k)) THEN
    vcheck = volfx(MineralId(k),jx,jy,jz) + dppt(k,jx,jy,jz)*volmol(MineralId(k))*delt
  ELSE
    vcheck = volfx(k,jx,jy,jz) + dppt(k,jx,jy,jz)*volmol(k)*delt
  END IF
  
  IF (vcheck < 0.0) THEN

    IF (MineralAssociate(k)) THEN
      dppt(k,jx,jy,jz) = -volfx(MineralID(k),jx,jy,jz)/(volmol(MineralId(k))*delt)
      rmin(1,k) = dppt(k,jx,jy,jz)
      IF (nreactmin(k) > 1) THEN
        DO np = 2,nreactmin(k)
          rmin(np,k) = 0.0
        END DO
      END IF
      ivolume(k) = 1

    ELSE

!!      IF (volfx(k,jx,jy,jz) == 0.d0 .AND. dppt(k,jx,jy,jz) < 0.0d0) THEN
!!        ivolume(k) = 1
!!        dppt(k,jx,jy,jz) = 0.0d0
!!      ELSE

        dppt(k,jx,jy,jz) = -volfx(k,jx,jy,jz)/(volmol(k)*delt)
        rmin(1,k) = dppt(k,jx,jy,jz)
        IF (nreactmin(k) > 1) THEN
          DO np = 2,nreactmin(k)
            rmin(np,k) = 0.0
          END DO
        END IF
        ivolume(k) = 1
!!      END IF

    END IF
    
  END IF
  
END DO     !  End of kth mineral

RETURN
END SUBROUTINE reaction
