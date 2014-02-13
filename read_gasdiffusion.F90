!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:44
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE read_gasdiffusion(nout,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz

!  Default values here, overridden below

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
  
!  Check "gas_diffusion" string which fixes diffusion coefficient (ignoring temperature)
  
  IF (ssch == 'gas_diffusion') THEN    ! Overrides other specifications of diffusion coefficient
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      dgas = 0.0
      RETURN
    END IF
    
!  Looking for diffusion coefficient
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        dgas = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "gas_diffusion" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "gas_diffusion" '
      WRITE(*,*) ' Assuming gas diffusion coefficient = 0'
      WRITE(*,*)
      dgas = 0.0
    END IF
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
   
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_gasdiffusion
