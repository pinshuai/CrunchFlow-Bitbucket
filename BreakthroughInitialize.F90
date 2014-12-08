subroutine BreakthroughInitialize(ncomp,nspec,nkin,nrct,ngas,npot, &
    nx,ny,nz,nseries,nexchange,nexch_sec,nsurf,nsurf_sec )

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
USE io
USE strings
USE ReadFlow
USE modflowModule
USE NanoCrystal

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nkin
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: ngas
INTEGER(I4B), INTENT(IN)                                     :: npot
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz
INTEGER(I4B), INTENT(IN)                                     :: nseries
INTEGER(I4B), INTENT(IN)                                      :: nexchange
INTEGER(I4B), INTENT(IN)                                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                                      :: nsurf
INTEGER(I4B), INTENT(IN)                                      :: nsurf_sec

!! INTERNAL VARIABLES

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: nex
INTEGER(I4B)                                                  :: nmin
INTEGER(I4B)                                                  :: ikph

INTEGER(I4B)                                                  :: j
INTEGER(I4B)                                                  :: nxyz
INTEGER(I4B)                                                  :: jx
INTEGER(I4B)                                                  :: jy
INTEGER(I4B)                                                  :: jz
INTEGER(I4B)                                                  :: ll
INTEGER(I4B)                                                  :: intfile

CHARACTER (LEN=mls)                                           :: filename
CHARACTER (LEN=mls)                                           :: dummy
CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=12)                                            :: writeph

LOGICAL(LGT)                                                  :: ext

INTEGER(I4B)                                                  :: nplotsurface
INTEGER(I4B)                                                  :: nplotexchange
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: ExchangeBasis
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: SurfaceBasis

writeph = 'pH           '

IF (nexchange > 0) THEN
  nplotexchange = ncomp
  ALLOCATE(ExchangeBasis(ncomp))
  DO i = 1,ncomp
    Dummy = 'EXCH-'
    dumstring = ulab(i)
    Dummy(6:mls) = dumstring(1:mls-5)
    ExchangeBasis(i) = dummy
  END DO
ELSE
  nplotexchange = 0
END IF

IF (nsurf > 0) THEN
  nplotsurface = ncomp
  ALLOCATE(SurfaceBasis(ncomp))
  DO i = 1,ncomp
    Dummy = 'SURF-'
    dumstring = ulab(i)
    Dummy(6:mls) = dumstring(1:mls-5)
    SurfaceBasis(i) = dummy
  END DO
ELSE
  nplotsurface = 0
END IF

IF (nseries >= 1) THEN
  DO ll = 1,nseries
    intfile = 100+ll  
    IF (irestart == 1 .AND. AppendRestart) THEN    !  Open breakthrough file and go to end of file
      filename = TimeSeriesFile(ll)
      INQUIRE(FILE=filename,EXIST=ext)
      IF (.NOT. ext) THEN
        CALL stringlen(filename,ls)
        WRITE(*,*) 
        WRITE(*,*) ' Cannot find time series file: ', filename(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702,POSITION='append')
    ELSE
      filename = TimeSeriesFile(ll)
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702) 
!!CIS: Moved to individual plotting options      WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
2283 FORMAT('# Time series at grid cell: 'i3,1x,i3,1x,i3)      

      IF (OutputTimeUnits == 'years') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3001) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3001) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3001) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3007) (TimeSeriesUnits(i),i=1,nplot)
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3701) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3701) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3701) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF
        
      ELSE IF (OutputTimeUnits == 'days') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)

          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3002) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3002) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF

          ELSE
            WRITE(intfile,3002) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF

        ELSE IF (originlab) THEN

          IF (JennyDruhan) THEN
            WRITE(intfile,3036) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3008) (TimeSeriesUnits(i),i=1,nplot) 
          ELSE
            WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3008) (TimeSeriesUnits(i),i=1,nplot) 
          END IF
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)

          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3702) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3702) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF    

          ELSE

!! Kaleidagraph Jenny Druhan
            if (JennyDruhan) then
!!              WRITE(intfile,3702) (TimeSeriesSpecies(i),i=1,nplot), 'delCaAqueous' , 'delCaMineralBulk', 'delCaMineralInstant'
              WRITE(intfile,3712) (TimeSeriesSpecies(i),i=1,nplot) 
            else
              WRITE(intfile,3702) (TimeSeriesSpecies(i),i=1,nplot) 
            end if

          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'hours') THEN

        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3003) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3003) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3003) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
            WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
            WRITE(intfile,3009) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3703) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3703) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3703) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'minutes') THEN
 
        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3004) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3004) (ulab(iplot(i)),i=1,nplot),                &
               (ulab(i),i=1,nplotsurface),(ulab(i),i=1,nplotexchange),       &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3004) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3010) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3704) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3704) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3704) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF

      ELSE IF (OutputTimeUnits == 'seconds') THEN
        IF (tecplot) THEN
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3005) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3005) (ulab(iplot(i)),i=1,nplot),                &
               (ulab(i),i=1,nplotsurface),(ulab(i),i=1,nplotexchange),       &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3005) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ELSE IF (originlab) THEN
          WRITE(intfile,3006) (TimeSeriesSpecies(i),i=1,nplot) 
          WRITE(intfile,3011) (TimeSeriesUnits(i),i=1,nplot)  
        ELSE
          WRITE(intfile,2283) jxseries(ll),jyseries(ll),jzseries(ll)
          IF (iplotall == 1) THEN
            IF (ikph /= 0) THEN
              WRITE(intfile,3705) writeph,(ulab(iplot(i)),i=1,nplot),        &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            ELSE
              WRITE(intfile,3705) (ulab(iplot(i)),i=1,nplot),                &
                (SurfaceBasis(i),i=1,nplotsurface),                          &
                (ExchangeBasis(i),i=1,nplotexchange),                        &
                ('a_'//ulab(ik),ik=1,ncomp+nspec),                                 &
                (namsurf(ns),ns=1,nsurf),(namsurf_sec(ns),ns=1,nsurf_sec),   &
                (nam_exchsec(nex),nex=1,nexch_sec),                          &
                (umin(nmin),nmin=1,nrct)
            END IF
          ELSE
            WRITE(intfile,3705) (TimeSeriesSpecies(i),i=1,nplot) 
          END IF
        ENDIF
      ELSE
        WRITE(*,*) 
        WRITE(*,*) ' Output time units not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

    END IF
  END DO
END IF

IF (ALLOCATED(SurfaceBasis)) THEN
  DEALLOCATE(SurfaceBasis)
END IF
IF (ALLOCATED(ExchangeBasis)) THEN
  DEALLOCATE(ExchangeBasis)
END IF

!! Tecplot formats
3001 FORMAT('VARIABLES = "Time (yrs) " ',                   100(', "',A19,'"'))
3002 FORMAT('VARIABLES = "Time (days)" ',                   100(', "',A19,'"'))
3003 FORMAT('VARIABLES = "Time (hrs) " ',                   100(', "',A19,'"'))
3004 FORMAT('VARIABLES = "Time (min) " ',                   100(', "',A19,'"'))
3005 FORMAT('VARIABLES = "Time (sec) " ',                   100(', "',A19,'"'))

3006 FORMAT('  Time          ',                                       100(A17) )
3007 FORMAT('  Yrs           ',                                        100(A17) )
3008 FORMAT('  Days          ',                                        100(A17) )
3009 FORMAT('  Hrs           ',                                        100(A17) )
3010 FORMAT('  Min           ',                                        100(A17) )
3011 FORMAT('  Sec           ',                                        100(A17) )

3036 FORMAT('  Time     ', 2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                          'delCaMineral        ',3x,     &
                                          'del34S_Sulfate      ',3x,     &
                                          'del34S_Sulfide      ',3x,     &
                                          'del34S_Mineral      '     )  
3712 FORMAT('  Time(day)',2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                          'delCaMineral        ',3x,     &
                                          'del34S_Sulfate      ',3x,     &
                                          'del34S_Sulfide      ',3x,     &
                                          'del34S_Mineral      '     )
!! Kaleidagraph formats
3701 FORMAT('  Time(yrs)',17x,150(1X,A22))
3702 FORMAT('  Time(day)',17x,150(1X,A22))
3703 FORMAT('  Time(hrs)',17x,150(1X,A22))
3704 FORMAT('  Time(min)',17x,150(1X,A22))
3705 FORMAT('  Time(sec)',17x,150(1X,A22))

RETURN

702 WRITE(*,*)
IF (ALLOCATED(SurfaceBasis)) THEN
  DEALLOCATE(SurfaceBasis)
END IF
IF (ALLOCATED(ExchangeBasis)) THEN
  DEALLOCATE(ExchangeBasis)
END IF
WRITE(*,*) ' Error opening breakthrough file'
WRITE(*,*)
READ(*,*)
STOP


END subroutine BreakthroughInitialize