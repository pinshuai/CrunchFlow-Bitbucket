subroutine AqueousFluxInitialize(ncomp,nspec,nx,ny,nz )

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
INTEGER(I4B), INTENT(IN)                                     :: nx
INTEGER(I4B), INTENT(IN)                                     :: ny
INTEGER(I4B), INTENT(IN)                                     :: nz

!! INTERNAL VARIABLES

INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: nex
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

IF (nAqueousFluxSeriesFile >= 1) THEN
  DO ll = 1,nAqueousFluxSeriesFile
    intfile = 50+ll  
    IF (irestart == 1 .AND. AppendRestart) THEN    !  Open breakthrough file and go to end of file
      filename = AqueousFluxSeriesFile(ll)
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

      filename = AqueousFluxSeriesFile(ll)
      OPEN(UNIT=intfile,FILE=filename,STATUS='unknown',ERR=702) 

2283 FORMAT('# Aqueous flux series at grid cells: 'i3,'-',i3,1x,i3,'-',i3,1x,i3,'-',i3)   
2284 FORMAT('# Flux in cumulative moles' )  

        IF (tecplot) THEN

          WRITE(intfile,2283) jxAqueousFluxSeries_lo(ll),jxAqueousFluxSeries_hi(ll), jyAqueousFluxSeries_lo(ll),jyAqueousFluxSeries_hi(ll),   &
                              jzAqueousFluxSeries_lo(ll),jzAqueousFluxSeries_hi(ll)
          WRITE(intfile,2284)
!!          WRITE(intfile,3001) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux), (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)
          WRITE(intfile,3001) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)
        ELSE IF (originlab) THEN

!!          WRITE(intfile,3006) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux), (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)
          WRITE(intfile,3006) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)
!!          WRITE(intfile,3007) (TimeSeriesUnits(i),i=1,nplotAqueousFlux)

        ELSE

          WRITE(intfile,2283) jxAqueousFluxSeries_lo(ll),jxAqueousFluxSeries_hi(ll), jyAqueousFluxSeries_lo(ll),jyAqueousFluxSeries_hi(ll),   &
                              jzAqueousFluxSeries_lo(ll),jzAqueousFluxSeries_hi(ll)
          WRITE(intfile,2284)
!!          WRITE(intfile,3701) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux), (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)
          WRITE(intfile,3701) (AqueousFluxSeriesSpecies(i),i=1,nplotAqueousFlux)

        ENDIF

    END IF
  END DO
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

!! Kaleidagraph formats
3701 FORMAT('  Time(yrs)',17x,150(1X,A22))
3702 FORMAT('  Time(day)',17x,150(1X,A22))
3703 FORMAT('  Time(hrs)',17x,150(1X,A22))
3704 FORMAT('  Time(min)',17x,150(1X,A22))
3705 FORMAT('  Time(sec)',17x,150(1X,A22))

RETURN

702 WRITE(*,*)
WRITE(*,*) ' Error opening Aqueous Flux file'
WRITE(*,*)
READ(*,*)
STOP


END subroutine AqueousFluxInitialize