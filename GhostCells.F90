SUBROUTINE GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,xxx,TEXT)
USE crunchtype

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: lowX
INTEGER(I4B), INTENT(IN)                                    :: lowY
INTEGER(I4B), INTENT(IN)                                    :: lowZ
INTEGER(I4B), INTENT(IN)                                    :: highX
INTEGER(I4B), INTENT(IN)                                    :: highY
INTEGER(I4B), INTENT(IN)                                    :: highZ

REAL(DP), DIMENSION(lowX:highX,lowY:highY,lowZ:highZ), INTENT(INOUT)                   :: xxx

CHARACTER (LEN=15), INTENT(IN)                              :: text

!  Internal variables and arrays

INTEGER(I4B)                                                :: nxcheck
INTEGER(I4B)                                                :: nycheck
INTEGER(I4B)                                                :: nzcheck
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                     :: xdummy

IF (lowX /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in X direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (lowY /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in Y direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (lowZ /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in Z direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highX /= nx+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in X direction should be nx+2 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highY /= ny+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in Y direction should be ny+2 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highZ /= nz+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in X direction should be nz+2 for: ', text
  WRITE(*,*) 
  STOP
END IF

xxx(-1,:,:) = xxx(1,:,:)
xxx(0,:,:) = xxx(1,:,:)
xxx(nx+2,:,:) = xxx(nx,:,:)
xxx(nx+1,:,:) = xxx(nx,:,:)

xxx(:,-1,:) = xxx(:,1,:)
xxx(:,0,:) = xxx(:,1,:)
xxx(:,ny+2,:) = xxx(:,ny,:)
xxx(:,ny+1,:) = xxx(:,ny,:)

xxx(:,:,-1) = xxx(:,:,1)
xxx(:,:,0) = xxx(:,:,1)
xxx(:,:,nz+2) = xxx(:,:,nz)
xxx(:,:,nz+1) = xxx(:,:,nz)

RETURN

END SUBROUTINE GhostCells
