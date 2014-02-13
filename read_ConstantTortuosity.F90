SUBROUTINE read_ConstantTortuosity(nout,nx,ny,nz,constant_tortuosity,TortuosityConstant,TortuosityOption)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
LOGICAL(LGT), INTENT(OUT)                                   :: constant_tortuosity
REAL(DP), INTENT(OUT)                                       :: TortuosityConstant
CHARACTER (LEN=mls), INTENT(INOUT)                          :: TortuosityOption

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz
constant_tortuosity = .FALSE.
TortuosityConstant = 1.0d0        !! Set the tortuosity to 1.0 as a default
TortuosityOption = 'none'

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
! Check to see if initial substring is "constant_tortuosity"
  
  IF (ssch == 'constant_tortuosity') THEN
    IF (nxyz == 1) THEN    
      TortuosityConstant = 1.0d0
      constant_tortuosity = .TRUE.
      RETURN
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        TortuosityConstant = DNUM(ssch)
        constant_tortuosity = .TRUE.
      ELSE                !  An ascii string--so read it in to as an option to calculate tortuosity from the porosity (check later to see if it is a valid option)
        TortuosityOption = ssch
        TortuosityConstant = 1.0d0
        constant_tortuosity = .TRUE.              !!  Although not constant during time stepping, this will cause the code to skip over a file read or tortuosity by zone
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "constant_tortuosity" '
      WRITE(*,*) ' Value(s) must be given if constant_tortuosity specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_ConstantTortuosity
