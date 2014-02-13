! *************** Subroutine modread *********************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-09-15  Time: 21:39:08
 
! Purpose is to read unformatted .hff file created by modflow lkmt2.f
! Written by Steve Yabusaki and Randal Taira (3/2000)
! ***********************************************************************

SUBROUTINE modread(nx,ny,nz,activeflow)

USE crunchtype
USE medium
USE transport
USE flow
USE modflowModule


IMPLICIT NONE

INTERFACE
  SUBROUTINE readtype1(ndim,iii,jjj,kkk,xxx,num)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                               :: ndim
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: iii
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: jjj
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT)              :: kkk
  REAL(DP), DIMENSION(:), INTENT(OUT)                    :: xxx
  INTEGER(I4B), INTENT(OUT)                              :: num
  END SUBROUTINE readtype1
END INTERFACE

!  ****************** END INTERFACE BLOCKS  ********************

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

LOGICAL, INTENT(INOUT)                                   :: activeflow

! Internal variables and arrays

REAL(DP), DIMENSION(:), ALLOCATABLE                      :: qdum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: idum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: jdum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: kdum

INTEGER(I4B)                                             :: kper
INTEGER(I4B)                                             :: kstp
INTEGER(I4B)                                             :: ncol
INTEGER(I4B)                                             :: nrow
INTEGER(I4B)                                             :: nlay
INTEGER(I4B)                                             :: kperold
INTEGER(I4B)                                             :: kstpold
INTEGER(I4B)                                             :: i
INTEGER(I4B)                                             :: ncnhcheck
INTEGER(I4B)                                             :: nwellscheck
INTEGER(I4B)                                             :: nriverscheck
INTEGER(I4B)                                             :: ndrainscheck
INTEGER(I4B)                                             :: ndumcheck
INTEGER(I4B), PARAMETER                                  :: ndimen = 300           

CHARACTER (LEN=16)                                       :: text

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: jxLoc
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: jyLoc
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                  :: jzLoc

ALLOCATE(idum(ndimen))
ALLOCATE(jdum(ndimen))
ALLOCATE(kdum(ndimen))
ALLOCATE(qdum(ndimen))
ALLOCATE(jxLoc(ndimen))
ALLOCATE(jyLoc(ndimen))
ALLOCATE(jzLoc(ndimen))

! Start reading thickness, X,Y,Z flows, constant head,
!   and source/sink information

READ(1)

!  kper is presumably the stress period, kstp is the timestep within the stress period

READ(1,END = 1000) kper,kstp,ncol,nrow,nlay,text
kperold=kper
kstpold=kstp

DO
  select CASE (text)

!   read thickness
    CASE ("THKSAT")
      CALL readtype3(nx,ny,nz,text,dzz(1:nx,1:ny,1:nz))
!   read X flows
    CASE ("QXX")
      CALL readtype3(nx,ny,nz,text,qx(1:nx,1:ny,1:nz))
!   read Y flows
    CASE ("QYY")
      CALL readtype3(nx,ny,nz,text,qy(1:nx,1:ny,1:nz))
!   read Z flows
    CASE ("QZZ")
      CALL readtype3(nx,ny,nz,text,qz(1:nx,1:ny,1:nz))

!   read constant head locations and rates
    CASE ("CNH")
      CALL readtype1(ncnh,jxLoc,jyLoc,jzLoc,qcnh,ncnhcheck)
      IF (ncnhcheck /= ncnh) THEN
        WRITE(*,*) 
        WRITE(*,*) ' Mismatch in dimensions of constant head source terms'
        WRITE(*,*) ' Dimension:                                      ', ncnh
        WRITE(*,*) ' Value returned from read of MODFLOW "hff" file: ', ncnhcheck
        WRITE(*,*) 
        READ(*,*)
        STOP
      END IF

      DO i = 1,ncnh
        IF (jxLoc(i) /= jxHeadLoc(i) .OR. jyLoc(i) /= jyHeadLoc(i)   &
            .OR. jzLoc(i) /= jzHeadLoc(i)) THEN
          WRITE(*,*) 
          WRITE(*,*)   ' Mismatch between locations of constant head in MODFLOW and Crunch input files'
          WRITE(*,*)   ' Constant head number: ',i
          WRITE(*,501) ' MODFLOW location: ',jxLoc(i),jyLoc(i),jzLoc(i)
          WRITE(*,502) ' CRUNCH location : ',jxHeadLoc(i),jyHeadLoc(i),jzHeadLoc(i)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END DO

!   read well locations and rates
    CASE ("WEL")
      CALL readtype1(nwells,jxLoc,jyLoc,jzLoc,q,nwellscheck)

      IF (nwellscheck /= nwells) THEN
        WRITE(*,*) 
        WRITE(*,*) ' Mismatch in dimensions of wells'
        WRITE(*,*) ' Dimension:                                      ', nwells
        WRITE(*,*) ' Value returned from read of MODFLOW "hff" file: ', nwellscheck
        WRITE(*,*) 
        READ(*,*)
        STOP
      END IF

      DO i = 1,nwells
        IF (jxLoc(i) /= jxWellLoc(i) .OR. jyLoc(i) /= jyWellLoc(i)   &
            .OR. jzLoc(i) /= jzWellLoc(i)) THEN
          WRITE(*,*) 
          WRITE(*,*)   ' Mismatch between locations of well in MODFLOW and Crunch input files'
          WRITE(*,*)   ' Well number: ',i
          WRITE(*,501) ' MODFLOW location: ',jxLoc(i),jyLoc(i),jzLoc(i)
          WRITE(*,502) ' CRUNCH location : ',jxWellLoc(i),jyWellLoc(i),jzWellLoc(i)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END DO

      501 FORMAT(1x,'MODFLOW location: ',i4,1x,i4,1x,i4)
      502 FORMAT(1x,'CRUNCH location : ',i4,1x,i4,1x,i4)

!   read drain locations and rates

    CASE ("DRN")
      CALL readtype1(ndrains,jxLoc,jyLoc,jzLoc,qdrain,ndrainscheck)

      IF (ndrainscheck /= ndrains) THEN
        WRITE(*,*) 
        WRITE(*,*) ' Mismatch in dimensions of drains'
        WRITE(*,*) ' Dimension:                                      ', ndrains
        WRITE(*,*) ' Value returned from read of MODFLOW "hff" file: ', ndrainscheck
        WRITE(*,*) 
        READ(*,*)
        STOP
      END IF

      DO i = 1,ndrains
        IF (jxLoc(i) /= jxDrainLoc(i) .OR. jyLoc(i) /= jyDrainLoc(i)   &
            .OR. jzLoc(i) /= jzDrainLoc(i)) THEN
          WRITE(*,*) 
          WRITE(*,*)   ' Mismatch between locations of drain in MODFLOW and Crunch input files'
          WRITE(*,*)   ' Drain number: ',i
          WRITE(*,501) ' MODFLOW location: ',jxLoc(i),jyLoc(i),jzLoc(i)
          WRITE(*,502) ' CRUNCH location : ',jxDrainLoc(i),jyDrainLoc(i),jzDrainLoc(i)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END DO

!   read river locations and rates
    CASE ("RIV")
      CALL readtype1(nrivers,jxLoc,jyLoc,jzLoc,qriver,nriverscheck)

      IF (nriverscheck /= nrivers) THEN
        WRITE(*,*) 
        WRITE(*,*) ' Mismatch in dimensions of drains'
        WRITE(*,*) ' Dimension:                                      ', nrivers
        WRITE(*,*) ' Value returned from read of MODFLOW "hff" file: ', nriverscheck
        WRITE(*,*) 
        READ(*,*)
        STOP
      END IF

      DO i = 1,nrivers
        IF (jxLoc(i) /= jxRiverLoc(i) .OR. jyLoc(i) /= jyRiverLoc(i)   &
            .OR. jzLoc(i) /= jzRiverLoc(i)) THEN
          WRITE(*,*) 
          WRITE(*,*)   ' Mismatch between locations of river in MODFLOW and Crunch input files'
          WRITE(*,*)   ' River number: ',i
          WRITE(*,501) ' MODFLOW location: ',jxLoc(i),jyLoc(i),jzLoc(i)
          WRITE(*,502) ' CRUNCH location : ',jxRiverLoc(i),jyRiverLoc(i),jzRiverLoc(i)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END DO

!   read recharge locations and rates
    CASE ("RCH")
      CALL readtype2(nx,ny,lrecharge,qrecharge)
!   read evapotranspiration locations and rates
    CASE ("EVT")
      CALL readtype2(nx,ny,levt,qevt)
!   read head dependent boundary and rates
    CASE ("GHB")
      CALL readtype1(ndimen,idum,jdum,kdum,qdum,ndumcheck)
    CASE default
!     GO TO 1000
  END select
  READ(1,END = 1000) kper,kstp,ncol,nrow,nlay,text

!   If timestep or stress period changes, backspace and exit read.  
!     Data will be read during next timestep

!  NOTE:  Stress period or timestep should not change as multiple "cases" are read from a
!         single time step (e.g., velocity, well information)

  IF(kper /= kperold.OR.kstp /= kstpold)THEN

!   Update ghost cell dzz

    dzz(0,:,:)=dzz(1,:,:)
    dzz(:,0,:)=dzz(:,1,:)
    dzz(:,:,0)=dzz(:,:,1)
    dzz(:,:,-1)=dzz(:,:,1)
    dzz(:,:,nz+1)=dzz(:,:,nz)
    dzz(:,:,nz+2)=dzz(:,:,nz)

!   Calculate common factors

!!    bydzz = 1.0d0 / dzz(1:nx,1:ny,1:nz)
!!    DO jz = 1, nz-1
!!      DO jy = 1, ny
!!        DO jx = 1, nx
!!          bydistz(jx,jy,jz) = 2.0d0/(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
!!        END DO
!!      END DO
!!    END DO

!   Classify source/sink terms

    DO i = 1, nwells
      IF (Q(i) < 0.0) THEN
        Wtype(i) = 'Extraction'
      ELSE
        Wtype(i) = 'Injection'
      END IF
    END DO

    DO i = 1, ncnh
      IF (Qcnh(i) < 0.0) THEN
        Htype(i) = 'Extraction'
      ELSE
        Htype(i) = 'Injection'
      END IF
    END DO

    DO i=1,nrivers
      IF (Qriver(i) < 0.0) THEN
        Rtype(i) = 'Extraction'
      ELSE
        Rtype(i) = 'Injection'
      END IF
    END DO

    BACKSPACE 1
    BACKSPACE 1
    DEALLOCATE(idum)
    DEALLOCATE(jdum)
    DEALLOCATE(kdum)
    DEALLOCATE(qdum)
    DEALLOCATE(jxLoc)
    DEALLOCATE(jyLoc)
    DEALLOCATE(jzLoc)
    RETURN
  END IF

END DO
1000 activeflow=.FALSE.
DEALLOCATE(idum)
DEALLOCATE(jdum)
DEALLOCATE(kdum)
DEALLOCATE(qdum)
DEALLOCATE(jxLoc)
DEALLOCATE(jyLoc)
DEALLOCATE(jzLoc)
RETURN
END SUBROUTINE modread

! **********  End Main Program modread  ********************************


!      SUBROUTINE getfile()
! ***********************************************************************
!  get hff filename from user
! ***********************************************************************

!      integer gridDims
! character*16 hffname

!      print *, "Enter hff file name:"
! read *, hffname
! open(1,file= hffname,form='unformatted', status = 'old')
! READ(1)

! return
! end

! ************* End SUBROUTINE getfile **********************************


!      SUBROUTINE divergence(nx,ny,nz,flowX,flowY,flowZ)
! ***********************************************************************
!  calculate divergence for interior cells and print to file
! ***********************************************************************

!      dimension flowX(nx,ny,nz)
!      dimension flowY(nx,ny,nz)
!      dimension flowZ(nx,ny,nz)
! dimension diverge(nx,ny,nz)

!      do k=1,nlay
!   do j=1,nrow
!         do i=1,ncol
!            diverge(i,j,k) = 0.0
!          end do
!   end do
! end do

!      write(2,*) 'Divergence numbers'

!      do k=2,nlay-1
!   do j=2,nrow-1
!          do i=2,ncol-1
!            diverge(i,j,k) = (flowX(i-1,j,k-1) - flowX(i,j,k-1))
!     &  + (flowY(i,j-1,k-1) - flowY(i,j,k-1)) +
!     &          (flowZ(i,j,k-1) - flowZ(i,j,k))
!          end do
!   end do
! end do

!   do k=1,nlay
!   do j=1,nrow
!     write(2,'(12(1pe12.5))')(diverge(i,j,k),i=1,ncol)
!   end do
! end do

! return
! end

! ************* End SUBROUTINE divergence *******************************



