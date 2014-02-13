SUBROUTINE read_TortuosityFile(nout,nx,ny,nz,constant_tortuosity,TortuosityFile,lfile,FileFormatType)
USE crunchtype
USE params
USE concentration
USE transport, ONLY: tortuosity
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
LOGICAL(LGT), INTENT(INOUT)                                 :: constant_tortuosity
CHARACTER (LEN=mls), INTENT(OUT)                            :: TortuosityFile
INTEGER(I4B), INTENT(OUT)                                   :: lfile
CHARACTER (LEN=mls), INTENT(OUT)                            :: FileFormatType

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat

nxyz = nx*ny*nz
TortuosityFile = ' '

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
  IF (ssch == 'read_tortuosity' .OR. ssch=='read_tortuosityfile') THEN
    IF (nxyz == 1) THEN    
      Tortuosity = 1.0d0
      constant_tortuosity = .TRUE.
      RETURN
    END IF
    constant_tortuosity = .FALSE.
    
!  Looking for tortuosity file name 
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL stringtype(ssch,lzs,res)
      TortuosityFile = ssch
      lfile = ls

!!  *****************************************************************
!!    Now, check for a file format
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      lenformat = ls
      CALL majuscules(ssch,lenformat)
      IF (ls /= 0) THEN
        IF (ssch == 'singlecolumn') THEN
          FileFormatType = 'SingleColumn'
        ELSE IF (ssch == 'continuousread') THEN
          FileFormatType = 'ContinuousRead'
        ELSE IF (ssch == 'unformatted') THEN
          FileFormatType = 'Unformatted' 
        ELSE IF (ssch == 'distanceplusvariable' .OR. ssch == 'fullform' .OR. ssch == 'full') THEN
          FileFormatType = 'FullForm'         
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' File format not recognized: ', ssch(1:lenformat)
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIf
      ELSE    !! No file format provided, so assume default
        FileFormatType = 'SingleColumn'
      ENDIF
!!  *****************************************************************

    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name following "read_TortuosityFile" in TRANSPORT section '
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
END SUBROUTINE read_tortuosityfile
