SUBROUTINE read_toperatio(nout,ncomp,ndecay,idecayall,idpoint,kd,nchem,topefound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: ndecay
INTEGER(I4B), INTENT(IN)                                    :: idecayall
INTEGER(I4B), INTENT(IN)                                    :: idpoint
INTEGER(I4B), INTENT(IN OUT)                                :: kd
INTEGER(I4B), INTENT(IN)                                    :: nchem
LOGICAL(LGT), INTENT(IN OUT)                                :: topefound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dumstring
CHARACTER (LEN=mls)                                         :: toperatio

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: ipoint
INTEGER(I4B)                                                :: nisot
INTEGER(I4B)                                                :: lcond
INTEGER(I4B)                                                :: ltope
INTEGER(I4B)                                                :: isotope

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: eps


eps = 1.e-10

ipoint = idecay(idpoint)
nisot = nisotope(idpoint)

dumstring = condlabel(nchem)
CALL stringlen(dumstring,lcond)
toperatio = ulab(ipoint)
CALL stringlen(toperatio,ltope)

allocate(namtope(nisot))
allocate(temp_ratio(nisot))

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'isotope_ratio') THEN    !  "isotope_ratio" label found
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs = ls
      IF (ssch == toperatio) THEN   !  Decay element found
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          IF (idecayall == 1) THEN
            IF (ssch == 'all' .OR. ssch == 'All' .OR. ssch == 'ALL') THEN
              OPEN(UNIT=nscratch,STATUS='scratch')
              REWIND nscratch
              WRITE(nscratch,*) zone(ids+ls+1:mls)
              REWIND nscratch
              READ(nscratch,*,ERR=500) (namtope(isotope),temp_ratio(isotope),  &
                  isotope=1,nisot)
              CLOSE(UNIT=nscratch)
              
!  Check first to see that ratios add up to 1
              
              sum = 0.0
              DO isotope = 1,nisot
                sum = sum + temp_ratio(isotope)
              END DO
              IF (sum < (1.0-eps) .OR. sum > (1.0+eps)) THEN
                WRITE(*,*)
                WRITE(*,*) ' Isotope ratios dont sum to 1'
                WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                WRITE(*,*)
                STOP
              END IF
              
!  Check to see that names of isotopes match what was given in the "DECAY" block
!  Assumes (for now) that the order is the same, otherwise...
              
              DO isotope=1,nisot
                IF (namtope(isotope) /= decay_label(isotope,idpoint)) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                  WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                  WRITE(*,*)
                  STOP
                END IF
              END DO
              topefound = .true.
              DO kd = 1,nmindecay(idpoint)
                DO isotope=1,nisot
                  ratio_isotope_init(isotope,kd,idpoint,nchem) = temp_ratio(isotope)
                END DO
              END DO
              
            ELSE     !  Label is not "all"
              GO TO 100
            END IF
            
          ELSE    !  "idecayall" equals 0, so don't search for "all" string
            
!  Check that we have the right mineral name
            
            IF (ssch == umin(kdecay(kd,idpoint))) THEN
              
              OPEN(UNIT=nscratch,STATUS='scratch')
              REWIND nscratch
              WRITE(nscratch,*) zone(ids+ls+1:mls)
              REWIND nscratch
              READ(nscratch,*,ERR=500) (namtope(isotope),temp_ratio(isotope),  &
                  isotope=1,nisot)
              CLOSE(UNIT=nscratch)
              
              sum = 0.0
              DO isotope = 1,nisot
                sum = sum + temp_ratio(isotope)
              END DO
              IF (sum < (1.0-eps) .OR. sum > (1.0+eps)) THEN
                WRITE(*,*)
                WRITE(*,*) ' Isotope ratios dont sum to 1'
                WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                WRITE(*,*)
                STOP
              END IF
              
!  Check to see that names of isotopes match what was given in the "DECAY" block
!  Assumes (for now) that the order is the same, otherwise...
              
              DO isotope=1,nisot
                IF (namtope(isotope) /= decay_label(isotope,idpoint)) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                  WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                  WRITE(*,*)
                  STOP
                END IF
              END DO
              topefound = .true.
              DO isotope=1,nisot
                ratio_isotope_init(isotope,kd,idpoint,nchem) = temp_ratio(isotope)
              END DO
            ELSE
              GO TO 100
            END IF
            
          END IF
          
          deallocate(namtope)
          deallocate(temp_ratio)
          RETURN    !  Information found, so exit subroutine
          
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Cannot find isotopic ratio for: ',toperatio(1:ltope)
          WRITE(*,*) ' In condition: ',dumstring(1:lcond)
          WRITE(*,*)
          STOP
        END IF
      END IF
      GO TO 100     !  Continue looking
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Label "isotope_ratio" should be followed by radioactive element name'
      WRITE(*,*) ' In condition: ',dumstring(1:lcond)
      WRITE(*,*)
      STOP
    END IF
    
    
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF


300 deallocate(namtope)
deallocate(temp_ratio)
RETURN

500 WRITE(*,*)
WRITE(*,*) ' Error reading isotopic ratio for: ',toperatio(1:ltope)
WRITE(*,*) ' In condition: ',dumstring(1:lcond)
WRITE(*,*)
STOP

END SUBROUTINE read_toperatio
